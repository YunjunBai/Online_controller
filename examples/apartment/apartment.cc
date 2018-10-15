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
const double tau = 0.3;

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
  state_type s_ub={{10,6,3.5}};
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

  double H[8][4] = {
    { 4  , 4.2, 0  ,   2 },
    { 5, 7.6, 1.8  , 2 },
    { 7.4, 7.6, 0  ,   2 },
    { 8.4, 10, 1.8  ,   2 },
    { 4, 4.2, 3.4  ,  6 },
    { 5, 7.6  , 3.4  , 3.6},
    { 7.4, 7.6  , 3.4  , 6 },
    { 8.4  , 10, 3.4  ,  3.6 }
  };
  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [&H,ss,s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    double c1= s_eta[0]/2.0+1e-10;
    double c2= s_eta[1]/2.0+1e-10;
    for(size_t i=0; i<8; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0] <= (H[i][1]+c1) && 
          (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };


 /* write obstacles to file */
  write_to_file(ss,avoid,"obstacles");

  disturbance_type w_1={{0.05, 0.1, 0.05}};
  disturbance_type w_2={{0.03, 0.1, 0.05}};
  disturbance_type w2_lb={{7.6-i_eta[0]*p1,0,-3.5}};
  disturbance_type w2_ub={{10,1.8+i_eta[1]*p1*21/38,3.5}};
  disturbance_type w_3={{0.03, 0.1, 0.05}};
  disturbance_type w3_lb={{4.2,0,-3.5}};
  disturbance_type w3_ub={{7.4,1.8,3.5}};
  double persent=(w2_ub[0]-w2_lb[0])*(w2_ub[1]-w2_lb[1])/((s_ub[0]-s_lb[0])*(s_ub[1]-s_lb[1]));
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
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    r_es=ge.gb_estimate(r,u);
    disturbance_type w = dis.get_disturbance(x,r_es,avoid);
   
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
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    r_es=ge.gb_estimate(r,u);
    disturbance_type w = dis.get_disturbance(x,r_es,avoid); 
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

  std::cout << "\nComputing the initial transition function (before distrubance changes): " << std::endl;
  /* transition function of symbolic model */
  scots::TransitionFunction tf_o;
  scots::TransitionFunction tf_o1d;
  scots::TransitionFunction tf_new;
  scots::TransitionFunction tf_standard;
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
  
  tt.tic();
  abs.compute_gb(tf_o,rs_post, avoid);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();
 std::cout << "Number of new transitions: " << tf_o.get_no_transitions() << std::endl;
  //if(!getrusage(RUSAGE_SELF, &usage))
 //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
  
  dis.update_disturbance(w_2, w2_lb, w2_ub,avoid);
  tt.tic();
  abs.compute_gb(tf_o1d,rs_post, avoid);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();
  std::cout << "Number of new transitions: " << tf_o1d.get_no_transitions() << std::endl;
  dis.update_disturbance(w_3, w3_lb, w3_ub,avoid);
  std::cout << "\nComputing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new,tf_o1d, w3_lb, w3_ub, rs_repost, avoid);
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  double t1=tt.toc();

   std::cout << "\nComputing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post,avoid);
  
 if(!getrusage(RUSAGE_SELF, &usage))
   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  double t2= tt.toc();
  
  

  /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (9 <= (x[0]-s_eta[0]/2.0) && (x[0]+s_eta[0]/2.0) <= 9.5 && 
        0 <= (x[1]-s_eta[1]/2.0) && (x[1]+s_eta[1]/2.0) <= 0.5)
      return true;
    return false;
  };
   /* write target to file */
  write_to_file(ss,target,"target");

 
  std::cout << "\nSynthesis: old controller" << std::endl;
  tt.tic();
  scots::WinningDomain win_1=scots::solve_reachability_game(tf_o1d,target);
  tt.toc();
  std::cout << "Winning domain size: " << win_1.get_size() << std::endl;
  
  std::cout << "\nWrite controller to controller_1.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win_1)),"controller_1"))
    std::cout << "Done. \n";

  std::cout << "\nSynthesis: new controller " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::solve_reachability_game(tf_new,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";


    std::ofstream write;
    std::ifstream read;
    write.open("result.txt", std::ios::app);          
    write << "p:"<<persent<<" t1:"<<t1<<" t2:"<<t2<<std::endl;
    write.close();
    read.close();

  // tt.tic();
  // std::queue<abs_type> online_queue = tf_o1d.get_difference(tf_new, win);
  // tt.toc();
  // std::cout<<"onlinequeuesize:"<<online_queue.size()<<std::endl;

  // std::cout << "\nOnline Synthesis: " << std::endl;
  // tt.tic();
  // scots::WinningDomain win_online=scots::online_reachability_game(tf_new, online_queue,avoid, win);
  // tt.toc();
  // std::cout << "Winning domain size: " << win_online.get_size() << std::endl;
  // std::cout << "\nWrite controller to online_controller.scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win_online)),"online_controller"))
  //   std::cout << "Done. \n";


 // std::cout<<"\nstatic synthesis:"<<std::endl;
 // scots::WinningDomain win_static=scots::static_reachability_game(tf_new,target);
 // std::cout << "\nWrite controller to controller_w2.scs \n";
 // if(write_to_file(scots::StaticController(ss,is,std::move(win_static)),"controller_w2"))
  //  std::cout << "Done. \n";
  
  
}
int  main()
{
 
  for (int k = 0; k < 39; ++k)
  {
    
      main_parameters(k);
  }
  return 1;
}
