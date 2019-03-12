/*
 * vehicle.cc
 *
 *  created: Oct 2016
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 * doi: 10.1109/TAC.2016.2593947
 */

#include <iostream>
#include <array>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"
#include "Simpson3.hh"
#include "MatrixExp.hh"


/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;

/* sampling time */
const double tau = 0.15;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using disturbance_type = std::array<double, state_dim>;
using input_type = std::array<double,input_dim>;
using ds_type = std::array<double, 2*state_dim>;
using matrix_type = std::array<std::array<double,state_dim>,state_dim>;
/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;


void main_parameters(const int p1){ 
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{0,0,-3.5}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{15,5,3.5}};
  /* grid node distance diameter */
  state_type s_eta={{.2,.2,.2}};
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{-1,-1}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{ 1, 1}};
  /* grid node distance diameter */
  input_type i_eta={{.3,.3}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  double H[5][4] = {
    { 2  ,3, 2  ,   3 },
    { 4.8, 5, 0  , 3.4 },
    { 7, 8, 3  ,  3.4 },
    { 9.8, 10, 3  ,   3.2 },
    { 10, 11.4, 3 , 3.2 }
   
  };
  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [&H,ss,s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    double c1= s_eta[0]/2.0+1e-10;
    double c2= s_eta[1]/2.0+1e-10;
    for(size_t i=0; i<5; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0] <= (H[i][1]+c1) && 
          (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };


 /* write obstacles to file */
  write_to_file(ss,avoid,"obstacles");
  //disturbance_type w_1={{0.1, 0.1, 0.1}};
  disturbance_type w_1={{0.9, 0.9,0.9}};
  disturbance_type w_2={{0.01, 0.01, 0.01}};
  disturbance_type w2_lb={{5,0,-3.5}};
  disturbance_type w2_ub={{10,5,3.5}};
  disturbance_type w_3={{0.02, 0.02, 0.02}};
  disturbance_type w3_lb={{10,0,-3.5}};
  disturbance_type w3_ub={{15,2,3.5}};
 
  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);

  auto l_matrix=[&is](const abs_type& input_id){
    matrix_type l;
    input_type u;
    is.itox(input_id,u);
    for (int i = 0; i < state_dim; ++i)
      for (int j = 0; j < state_dim; ++j){
        if ((i==0&&j==2) ||(i==1&&j==2))
          l[i][j]=std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
        else
          l[i][j]=0;  
      }
    return l;
  };

  scots::GbEstimation<disturbance_type,matrix_type> ge(is, ss,w_1,w_2);
  ge.exp_interals(l_matrix,tau/10);
  
  auto rs_post = [&dis,&ge,avoid,w3_ub,w3_lb](ds_type &y, input_type &u) -> void {
   // dis.set_out_of_domain();
 
  auto rhs =[&dis,&ge,avoid](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    state_type r_es;
    disturbance_type w;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    
    r_es=ge.gb_estimate(r,u);
    w = dis.get_disturbance(x,r_es,avoid);
  

    //disturbance_type w = {0.05, 0.05, 0.05};
    double alpha=std::atan(std::tan(u[1])/2.0);
    double c = std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
    yy[0] = u[0]*std::cos(alpha+y[2])/std::cos(alpha);
    yy[1] = u[0]*std::sin(alpha+y[2])/std::cos(alpha);
    yy[2] = u[0]*std::tan(u[1]);
    yy[3] = c*y[5] + w[0];
    yy[4] = c*y[5] + w[1];
    yy[5] = 0;
  };
  //while(ignore==false){
    scots::runge_kutta_fixed4(rhs,y,u,dis, w3_lb,w3_ub, 2*state_dim,tau,10);
  //  ignore = dis.get_out_of_domain();
};

auto rs_repost = [&dis,&ge,w3_lb,w3_ub,avoid](ds_type &y, input_type &u, bool &neigbour) -> void {
  dis.set_intersection_check();
  
  //dis.set_out_of_domain();
  auto rhs =[&dis,&ge,w3_lb,w3_ub,avoid](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    state_type r_es;
    disturbance_type w;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    
      r_es=ge.gb_estimate(r,u);
    w = dis.get_disturbance(x,r_es,avoid); 
    
    
    
    //disturbance_type w = {0.05, 0.05, 0.05};
    double alpha=std::atan(std::tan(u[1])/2.0);
    double c = std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
    yy[0] = u[0]*std::cos(alpha+y[2])/std::cos(alpha);
    yy[1] = u[0]*std::sin(alpha+y[2])/std::cos(alpha);
    yy[2] = u[0]*std::tan(u[1]);
    yy[3] = c*y[5] + w[0];
    yy[4] = c*y[5] + w[1];
    yy[5] = 0;
   
  };
  //ignore = dis.get_out_of_domain();
  //if(ignore==false)    
  scots::runge_kutta_fixed4(rhs,y,u,dis,w3_lb,w3_ub, 2*state_dim,tau,10);

  if(dis.get_intersection_check()==true){
    neigbour=true;
  }
  
};

 /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (13 <= (x[0]-s_eta[0]/2.0) && (x[0]+s_eta[0]/2.0) <= 15 && 
        2 <= (x[1]-s_eta[1]/2.0) && (x[1]+s_eta[1]/2.0) <= 4)
      return true;
    return false;
  };
   /* write target to file */
  write_to_file(ss,target,"target");

  std::cout << "\nComputing the initial transition function (before distrubance changes): " << std::endl;
  /* transition function of symbolic model */
  scots::TransitionFunction tf_1;
  scots::TransitionFunction tf_2;
  scots::TransitionFunction tf_3;
  scots::TransitionFunction tf_4;
  scots::Abstraction<state_dim,input_dim> abs(ss,is);
  //std::queue<abs_type> online_queue; 

  
  tt.tic();
  abs.compute_gb(tf_1,rs_post, avoid);
  tt.toc();
 std::cout << "Number of new transitions: " << tf_1.get_no_transitions() << std::endl;
  if(!getrusage(RUSAGE_SELF, &usage))
   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_1.get_no_transitions() << std::endl;
  
    std::cout << "\nSynthesis:  controller_1" << std::endl;
  tt.tic();
  scots::WinningDomain win_1=scots::solve_reachability_game(tf_1,target,avoid);
  tt.toc();
  std::cout << "Winning domain size: " << win_1.get_size() << std::endl;
  std::cout << "\nWrite controller to controller_1.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win_1)),"controller_global_1"))
    std::cout << "Done. \n";


  dis.update_disturbance(w_2, w2_lb, w2_ub,avoid);
  
  std::cout << "\nComputing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_2, tf_1, false,w2_lb, w2_ub, rs_repost, avoid);
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_2.get_no_transitions() << std::endl;
  std::cout << "Number of new transitions: " << tf_2.get_no_transitions() << std::endl;
  tt.toc();
 // std::cout<<"onlinequeuesize:"<<online_queue.size()<<std::endl;

    std::cout << "\nSynthesis:  controller_2" << std::endl;
  tt.tic();
  scots::WinningDomain win_2=scots::solve_reachability_game(tf_2,target,avoid);
  tt.toc();
  std::cout << "Winning domain size: " << win_2.get_size() << std::endl;
  std::cout << "\nWrite controller to controller_2.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win_2)),"controller_global_2"))
    std::cout << "Done. \n";

  dis.update_disturbance(w_3, w3_lb, w3_ub,avoid);
  
  std::cout << "\nComputing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_3, tf_2, true,w3_lb, w3_ub, rs_repost, avoid);
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_3.get_no_transitions() << std::endl;
  std::cout << "Number of new transitions: " << tf_3.get_no_transitions() << std::endl;
  tt.toc();
 // std::cout<<"onlinequeuesize:"<<online_queue.size()<<std::endl;

    std::cout << "\nSynthesis:  controller_3" << std::endl;
  tt.tic();
  scots::WinningDomain win_3=scots::solve_reachability_game(tf_3,target,avoid);
  tt.toc();
  std::cout << "Winning domain size: " << win_3.get_size() << std::endl;
  std::cout << "\nWrite controller to controller_3.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win_3)),"controller_global_3"))
    std::cout << "Done. \n";
  
  

 

  // std::cout << "\nSynthesis: old controller" << std::endl;
  // tt.tic();
  // scots::WinningDomain win_o=scots::solve_reachability_game(tf_o,target);
  // tt.toc();
  // std::cout << "Winning domain size: " << win_o.get_size() << std::endl;
  // std::cout << "\nWrite controller to controller_1.scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win_o)),"controller_scots_2"))
  //   std::cout << "Done. \n";

  
  
  // std::cout << "\nWrite controller to controller.scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller_local"))
  //   std::cout << "Done. \n";


    // std::ofstream write;
    // std::ifstream read;
    // write.open("result.txt", std::ios::app);          
    // write << "p:"<<persent<<" t1:"<<t1<<" t2:"<<t2<<std::endl;
    // write.close();
    // read.close();



 // std::cout<<"\nstatic synthesis:"<<std::endl;
 // scots::WinningDomain win_static=scots::static_reachability_game(tf_new,target);
 // std::cout << "\nWrite controller to controller_w2.scs \n";
 // if(write_to_file(scots::StaticController(ss,is,std::move(win_static)),"controller_w2"))
  //  std::cout << "Done. \n";
  
  
}
int  main()
{
  /*
  for (int k = 0; k < 39; ++k)
  {
      main_parameters(k);
  }
  */
  main_parameters(20);
  return 0;
}
