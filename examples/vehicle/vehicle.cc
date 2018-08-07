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

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* we integrate the vehicle ode by tau sec (the result is stored in x)  */
// auto  vehicle_post = [](state_type &x, const input_type &u) {
//   /* the ode describing the vehicle */
//   auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) {
//     double alpha=std::atan(std::tan(u[1])/2.0);
//     xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
//     xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
//     xx[2] = u[0]*std::tan(u[1]);
//   };
//   /* simulate (use 10 intermediate steps in the ode solver) */
//   scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
// };

/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
// auto radius_post = [](state_type &r, const state_type &, const input_type &u, const disturbance_type &w) {
//   double c = std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
//   r[0] = r[0]+c*r[2]*tau + w[0]*tau;
//   r[1] = r[1]+c*r[2]*tau + w[1]*tau;
// };

int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{0,0,-3.5}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{10,10,3.5}};
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

  disturbance_type w_1={{0.05, 0.05, 0.05}};
  disturbance_type w_2={{0.03, 0.1, 0.08}};
  disturbance_type w2_lb={{3,3,0.5}};
  disturbance_type w2_ub={{5,5,1.5}};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);


  /* set up constraint functions with obtacles */
  double H[15][4] = {
    { 1  , 1.2, 0  ,   9 },
    { 2.2, 2.4, 0  ,   5 },
    { 2.2, 2.4, 6  ,  10 },
    { 3.4, 3.6, 0  ,   9 },
    { 4.6, 4.8, 1  ,  10 },
    { 5.8, 6  , 0  ,   6 },
    { 5.8, 6  , 7  ,  10 },
    { 7  , 7.2, 1  ,  10 },
    { 8.2, 8.4, 0  ,  8.5},
    { 8.4, 9.3, 8.3,  8.5},
    { 9.3, 10 , 7.1,  7.3},
    { 8.4, 9.3, 5.9,  6.1},
    { 9.3, 10 , 4.7,  4.9},
    { 8.4, 9.3, 3.5,  3.7},
    { 9.3, 10 , 2.3,  2.5}
  };

  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [&H,ss,s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    double c1= s_eta[0]/2.0+1e-10;
    double c2= s_eta[1]/2.0+1e-10;
    for(size_t i=0; i<15; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0] <= (H[i][1]+c1) && 
          (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };
 
  /* write obstacles to file */
  write_to_file(ss,avoid,"obstacles");
  

  auto rs_post = [&dis](ds_type &y, input_type &u) -> void {
   // dis.set_out_of_domain();
  auto rhs =[&dis](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    disturbance_type w = dis.get_disturbance(x,r);
     
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
    scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,tau,10);
  //  ignore = dis.get_out_of_domain();
  
   
};

auto rs_repost = [&dis,w2_lb,w2_ub](ds_type &y, input_type &u, bool &neigbour) -> void {
  dis.set_intersection_check();
  //dis.set_out_of_domain();
  auto rhs =[&dis,w2_lb,w2_ub](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    disturbance_type w = dis.get_disturbance(x,r);
   
    dis.intersection(x,r, w2_lb,w2_ub);
    
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
  scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,tau,10);

  if(dis.get_intersection_check()==true){
    neigbour=true;
  }
  
};



  std::cout << "\nComputing the initial transition function (before distrubance changes): " << std::endl;
  /* transition function of symbolic model */
  scots::TransitionFunction tf_o1d;
  scots::TransitionFunction tf_new;
  scots::TransitionFunction tf_standard;
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
  
  tt.tic();
  abs.compute_gb(tf_o1d,rs_post, avoid);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();

  //if(!getrusage(RUSAGE_SELF, &usage))
 //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
 // std::cout << "Number of transitions: " << tf_o1d.get_no_transitions() << std::endl;
  dis.update_disturbance(w_2, w2_lb, w2_ub);

   std::cout << "\nComputing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post,avoid);
  
 if(!getrusage(RUSAGE_SELF, &usage))
   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  tt.toc();
  
  std::cout << "\nComputing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new,tf_o1d,tf_standard, w2_lb, w2_ub, rs_repost, avoid);
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  tt.toc();

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

 
  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::static_reachability_game(tf_o1d,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";


 
 

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
  
  return 1;
}
