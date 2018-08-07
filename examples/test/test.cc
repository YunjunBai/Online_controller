/*
 * unicycle.cc
 *
 *  
 */
#include <iostream>
#include <fstream>

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


void main_parameters(const int p1, const int p2) {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type s_lb={0, 0, -M_PI-0.4};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={11.4, 11.4, M_PI+0.4};
  /* grid node distance diameter */
  state_type s_eta={0.6, 0.6, 0.3};
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={-1.2, -1.6};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={1.2, 1.6};
  /* grid node distance diameter */
  input_type i_eta={.3, .2};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  disturbance_type w_1={0.05, 0.05, 0.05};
  disturbance_type w_2={0.03, 0.1, 0.08};
  disturbance_type w2_lb={7-2*p1*s_eta[0], 5-2*p1*s_eta[1], -M_PI-0.4};
  disturbance_type w2_ub={8.2+2*p1*s_eta[0], 6.2+2*p1*s_eta[1], M_PI+0.4};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);



  /* avoid function returns 1 if x is in avoid set  */
  auto avoid = [](const abs_type& idx) {
    
    return false;
  };
 
  /* write obstacles to file */
  write_to_file(ss,avoid,"obstacles");
  

  auto rs_post = [&dis, p2](ds_type &y, input_type &u) -> void {
  auto rhs =[&dis, p2](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    disturbance_type w = dis.get_disturbance(x,r);

    yy[0] = u[0]*std::cos(y[2])*p2;
    yy[1] = u[0]*std::sin(y[2])*p2;
    yy[2] = u[1];
    yy[3] = y[5]*std::abs(u[0]) + w[0];
    yy[4] = y[5]*std::abs(u[0]) + w[1];
    yy[5] = 0;
  };

  scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,tau,10);
};

auto rs_repost = [&dis,w2_lb,w2_ub, p2](ds_type &y, input_type &u, bool &neigbour) -> void {
  dis.set_intersection_check();
  //dis.set_out_of_domain();
  auto rhs =[&dis,w2_lb,w2_ub, p2](ds_type &yy, const ds_type &y, input_type &u) -> void {
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
    yy[0] = u[0]*std::cos(y[2])*p2;
    yy[1] = u[0]*std::sin(y[2])*p2;
    yy[2] = u[1];
    yy[3] = y[5]*std::abs(u[0]) + w[0];
    yy[4] = y[5]*std::abs(u[0]) + w[1];
    yy[5] = 0;
   
  };
  //ignore = dis.get_out_of_domain();
  //if(ignore==false)    
  scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,tau,10);

  if(dis.get_intersection_check()==true){
    neigbour=true;
  }
  
};

  std::cout << "Computing the initial transition function (before distrubance changes): " << std::endl;
  /* transition function of symbolic model */
  scots::TransitionFunction tf_o1d,tf_new,tf_standard;
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
  
  tt.tic();
  abs.compute_gb(tf_o1d,rs_post, avoid);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  double t1=tt.toc();

  //if(!getrusage(RUSAGE_SELF, &usage))
 //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
 // std::cout << "Number of transitions: " << tf_o1d.get_no_transitions() << std::endl;
  dis.update_disturbance(w_2, w2_lb, w2_ub);

   std::cout << "Computing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post,avoid);
  
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  double t3=tt.toc();

  std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new,tf_o1d, w2_lb, w2_ub, rs_repost, avoid);
 
   std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  double t2=tt.toc();
  /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (10.79 <= (x[0]-s_eta[0]/2.0) && (x[0]+s_eta[0]/2.0) <= 11.3 && 
        0.1 <= (x[1]-s_eta[1]/2.0) && (x[1]+s_eta[1]/2.0) <= 0.61)
      return true;
    return false;
  };
   /* write target to file */
  write_to_file(ss,target,"target");

 
  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::static_reachability_game(tf_new,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";

  std::ofstream write;
  std::ifstream read;
  write.open("result.txt", std::ios::app);          
  write << "p1:"<<p1<<" p2:"<<p2<<" t1:"<<t1<<" t2:"<<t2<<" t3:"<<t3<<std::endl;
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
  // scots::WinningDomain win_static=scots::solve_reachability_game(tf_new,target);
  // std::cout << "\nWrite controller to controller_w2.scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win_static)),"controller_w2"))
  //   std::cout << "Done. \n";
  
}


int  main()
{
 
  for (int k = 0; k < 5; ++k)
  {
    for(int r = 1; r<6; r++)
      main_parameters(k,r);
  }
  return 1;
}

