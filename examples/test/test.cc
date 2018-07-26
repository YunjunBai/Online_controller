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
const int state_dim=2;
/* input space dim */
const int input_dim=2;

/* sampling time */
const double tau = 0.9;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;
using disturbance_type = std::array<double,state_dim>;
using ds_type = std::array<double, 2*state_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;


int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type s_lb={0,0};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={6,1.8};
  /* grid node distance diameter */
  state_type s_eta={0.6,0.6};
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={-1.3,-1.3};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={1.3, 1.3};
  /* grid node distance diameter */
  input_type i_eta={0.5, 0.5};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  

  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    if (2.1 <= (x[0]) && (x[0]) <= 2.3 &&
      0 <= (x[1]) && (x[1]) <= 1.2)
      return true;
    return false;
  };

  /* write obstacles to file */
  write_to_file(ss,avoid,"obstacles");

  disturbance_type w_1={0.05, 0.05};
  disturbance_type w_2={0.03, 0.1};
  disturbance_type w2_lb={0,0};
  disturbance_type w2_ub={5,1.8};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);


  auto rs_post = [&dis](ds_type &y, input_type &u, bool check_i, disturbance_type lb, disturbance_type ub) -> void {
  auto rhs =[&dis](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    disturbance_type w = dis.get_disturbance(x,r);
    dis.check_xinR(x,check_i,lb,ub);//todo
    double L[2][2];
        L[0][0] = 0;
        L[0][1] = 0;
        L[1][0] = 0;
        L[1][1] = 0;
    /* coupled system + growth bound ode */
    yy[0] = u[0];
    yy[1] = u[1];
    yy[2] = L[0][0]*y[2] + L[0][1]*y[2] + w[0];
    yy[3] = L[1][0]*y[2] + L[1][0]*y[3] + w[1];
  };
  scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,0.9,10);
};

  std::cout << "Computing the transition function: " << std::endl;
  /* transition function of symbolic model */
  scots::TransitionFunction tf_o1d,tf_new;
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
 
  tt.tic();
  abs.compute_gb(tf_o1d,rs_post, avoid);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  std::cout << "Time to compute transition: "<< "\n" <<std::endl;
  tt.toc();
  

  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_o1d.get_no_transitions() << std::endl;

  /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (2.5 <= (x[0]) && (x[0]) <= 3.5 &&
        -0.7 <= (x[1]) && (x[1]) <= 0.7)
      return true;
    return false;
  };
   /* write target to file */
  write_to_file(ss,target,"target");

 
  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::static_reachability_game(tf_o1d,target);
  tt.toc();
  std::cout << "Winning domain size (first step of online): " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  //if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
   // std::cout << "Done. \n";

  tt.tic();

  dis.update_disturbance(w_2, w2_lb, w2_ub);
  abs.compute_gb(tf_new,rs_post, avoid);
  std::cout << "Number of transitions: " << tf_new.get_no_transitions() << std::endl;
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();

 
  tt.tic();
  std::queue<abs_type> online_queue = tf_o1d.get_difference(tf_new, win);
  tt.toc();
   std::cout<<"onlinequeuesize:"<<online_queue.size()<<std::endl;

  std::cout << "\nOnline Synthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win_online=scots::online_reachability_game(tf_new, online_queue,avoid, win);
  tt.toc();
  std::cout << "Winning domain size (online algorithm): " << win_online.get_size() << std::endl;
  std::cout << "\nWrite controller to online_controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win_online)),"online_controller"))
    std::cout << "Done. \n";


  std::cout<<"\nstatic synthesis:"<<std::endl;
  tt.tic();
  scots::WinningDomain win_static=scots::solve_reachability_game(tf_new,target);
   tt.toc();
   std::cout << "Winning domain size (static algorithm): " << win_static.get_size() << std::endl;
  std::cout << "\nWrite controller to controller_w2.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win_static)),"controller_w2"))
    std::cout << "Done. \n";
  
  return 1;
}
