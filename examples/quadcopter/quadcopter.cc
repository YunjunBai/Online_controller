/*
 * quadcopter.cc
 *
 *  created: Sep 2018
 *  author: Yunjun
 *
 */

/*
 * information about this example is given in
 * http://charlestytler.com/quadcopter-equations-motion/
 *
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
const int state_dim=6;
/* input space dim */
const int input_dim=6;
/* sampling time */
const double tau = 0.25;
using abs_type = scots::abs_type;
/* data types of the state space elements and input 
 * space elements used in uniform grid and ode solver */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;
using disturbance_type = std::array<double, state_dim>;
using ds_type = std::array<double, 2*state_dim>;


int main() {
  /* to measure time */
  TicToc tt;

  /* construct grid for the state space */
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  /* optimized values computed according to doi: 10.1109/CDC.2015.7403185 */
  state_type s_eta={{1,1,1,20*M_PI/180,20*M_PI/180,20*M_PI/180}}; 
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{0,0,0,-20*M_PI/180,-20*M_PI/180,-20*M_PI/180}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{5,5,5,60*M_PI/180,60*M_PI/180,60*M_PI/180}}; 
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();

  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{0,0,0,0,0,0}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{0.5,0.5,0.5,8*M_PI/180,8*M_PI/180,8*M_PI/180}};
  /* grid node distance diameter */
  input_type i_eta={{.25,.25,.25, 2*M_PI/180, 2*M_PI/180, 2*M_PI/180}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* setup object to compute the transition function */
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
  /* measurement disturbances  */
  state_type z={{0.05,0.05,0.05,0.0025/180*M_PI,0.0025/180*M_PI,0.0025/180*M_PI}};
  abs.set_measurement_error_bound(z);

  disturbance_type w_1={.108,0.002,0,0,0,0};
  disturbance_type w_2={0.203, 0.001, 0.001,0,0,0};
  disturbance_type w2_lb={0,0,0,-20*M_PI/180,-20*M_PI/180,-20*M_PI/180};
  disturbance_type w2_ub={3,3,3,60*M_PI/180,60*M_PI/180,60*M_PI/180};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);

  auto rs_post = [&dis,w2_lb,w2_ub](ds_type &y, input_type &u) -> void {
  auto rhs_1 =[&dis](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    disturbance_type w = dis.get_disturbance(x,r);
    double L[6][6];
    L[0][3]=u[1]+u[2];
    L[0][4]=u[2]+0.9129*u[0]+u[1];
    L[0][5]=u[2]+0.5944*u[1]+0.1665*u[0];
    L[1][3]=u[1]+u[2];
    L[1][4]=u[2]+8.0165e-15*u[0]+0.8335*u[1];
    L[1][5]=1.0409*u[2]+0.9860*u[1]+u[0];
    L[2][3]=u[1]+0.9129*u[2];
    L[2][4]=u[0]+u[2]+u[1]*0.8335;
    L[3][3]=(3.3723e+11) * u[4] +(9.6712e+08)*u[5];
    L[3][4]=6.0049*u[5]+u[4]*(4.7596e+13);
    L[4][3]=u[4]+u[5];
    L[5][3]=2.4505*u[4]+2.2372*u[5];
    L[5][4]=u[4]*5.0049+u[5]*(8.5318+14);
    
    yy[0] = std::cos(y[4])*std::cos(y[5])*u[0]+(-std::cos(y[3])*std::sin(y[5])+std::sin(y[3])*std::sin(y[4])*std::cos(y[5]))*u[1]+(std::sin(y[3])*std::sin(y[5])+std::cos(y[3])*std::sin(y[4])*std::cos(y[5]))*u[3];
    yy[1] = std::cos(y[4])*std::sin(y[5])*u[0]+(std::cos(y[3])*std::cos(y[5])+std::sin(y[3])*std::sin(y[4])*std::sin(y[5]))*u[1]+(-std::sin(y[3])*std::cos(y[5])+std::cos(y[3])*std::sin(y[4])*std::sin(y[5]))*u[3];
    yy[2] = -std::sin(y[4])*u[0]+std::sin(y[3])*std::cos(y[4])*u[1]+std::cos(y[3])*std::cos(y[4])*u[2];
    yy[3] = u[3] +(u[4]*std::sin(y[3])+u[5]*std::cos(y[3]))*std::tan(y[4]);
    yy[4] = u[4]*std::cos(y[3])-u[5]*std::sin(y[3]);
    yy[5] = (u[4]*std::sin(y[3])+u[5]*std::cos(y[3]))*std::asin(y[4]);

    /* to account for input disturbances */
    yy[6] = L[0][3]*y[9]+L[0][4]*y[10]+L[0][5]*y[11]+w[0]; 
    yy[7] = L[1][3]*y[9]+L[1][4]*y[10]+L[1][5]*y[11]+w[1]; 
    yy[8] = L[2][3]*y[9]+L[2][4]*y[10]+w[2]; 
    yy[9] = L[3][3]*y[9]+L[3][4]*y[10]+w[3];
    yy[10]= L[4][3]*y[9]+w[4];
    yy[11]= L[5][3]*y[9]+L[5][4]*y[10]+w[5];
  };
  scots::runge_kutta_fixed4(rhs_1,y,u,dis, w2_lb,w2_ub,2*state_dim,tau,10);
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
   
    double L[6][6];
    L[0][3]=u[1]+u[2];
    L[0][4]=u[2]+0.9129*u[0]+u[1];
    L[0][5]=u[2]+0.5944*u[1]+0.1665*u[0];
    L[1][3]=u[1]+u[2];
    L[1][4]=u[2]+8.0165e-15*u[0]+0.8335*u[1];
    L[1][5]=1.0409*u[2]+0.9860*u[1]+u[0];
    L[2][3]=u[1]+0.9129*u[2];
    L[2][4]=u[0]+u[2]+u[1]*0.8335;
    L[3][3]=(3.3723e+11) * u[4] +(9.6712e+08)*u[5];
    L[3][4]=6.0049*u[5]+u[4]*(4.7596e+13);
    L[4][3]=u[4]+u[5];
    L[5][3]=2.4505*u[4]+2.2372*u[5];
    L[5][4]=u[4]*5.0049+u[5]*(8.5318+14);

    yy[0] = std::cos(y[4])*std::cos(y[5])*u[0]+(-std::cos(y[3])*std::sin(y[5])+std::sin(y[3])*std::sin(y[4])*std::cos(y[5]))*u[1]+(std::sin(y[3])*std::sin(y[5])+std::cos(y[3])*std::sin(y[4])*std::cos(y[5]))*u[3];
    yy[1] = std::cos(y[4])*std::sin(y[5])*u[0]+(std::cos(y[3])*std::cos(y[5])+std::sin(y[3])*std::sin(y[4])*std::sin(y[5]))*u[1]+(-std::sin(y[3])*std::cos(y[5])+std::cos(y[3])*std::sin(y[4])*std::sin(y[5]))*u[3];
    yy[2] = -std::sin(y[4])*u[0]+std::sin(y[3])*std::cos(y[4])*u[1]+std::cos(y[3])*std::cos(y[4])*u[2];
    yy[3] = u[3] +(u[4]*std::sin(y[3])+u[5]*std::cos(y[3]))*std::tan(y[4]);
    yy[4] = u[4]*std::cos(y[3])-u[5]*std::sin(y[3]);
    yy[5] = (u[4]*std::sin(y[3])+u[5]*std::cos(y[3]))*std::asin(y[4]);

    /* to account for input disturbances */
    yy[6] = L[0][3]*y[9]+L[0][4]*y[10]+L[0][5]*y[11]+w[0]; 
    yy[7] = L[1][3]*y[9]+L[1][4]*y[10]+L[1][5]*y[11]+w[1]; 
    yy[8] = L[2][3]*y[9]+L[2][4]*y[10]+w[2];
    yy[9] = L[3][3]*y[9]+L[3][4]*y[10]+w[3];
    yy[10]= L[4][3]*y[9]+w[4];
    yy[11]= L[5][3]*y[9]+L[5][4]*y[10]+w[5];
  };
  //ignore = dis.get_out_of_domain();
  //if(ignore==false)    
  scots::runge_kutta_fixed4(rhs,y,u,dis, w2_lb,w2_ub,2*state_dim,tau,10);

  if(dis.get_intersection_check()==true){
    neigbour=true;
  }
};

 /* transition function of symbolic model */
  scots::TransitionFunction tf_old, tf_new, tf_standard,tf_new_com;

  std::cout << "Computing the transition function: " << std::endl;
  tt.tic();
  abs.compute_gb(tf_old, rs_post);
  tt.toc();
  // if(!getrusage(RUSAGE_SELF, &usage))
  //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
  // std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

   dis.update_disturbance(w_2, w2_lb, w2_ub);
  //  double max_1=(1.0/60000) * (i_ub[0]-(2.7+3.08*std::pow(1.25+4.2*i_lb[1],2))*std::pow(s_lb[0],2)+ 60000.0*9.81);
  //  double max_2=(1/(60000*s_lb[0])) * (i_ub[0] + (68.6*(1.25+4.2*i_ub[1]))*std::pow(s_ub[0],2) + 60000.0*9.81);
  //  double max_3=s_ub[0];
  // state_type max_dynamic = {{ max_1, max_2 , max_3}};
  // state_type distance = dis.get_maxdistance(max_dynamic,tau);

   std::cout << "Computing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post);
  
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_standard.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  tt.toc();

  std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new,tf_old, w2_lb, w2_ub, rs_repost);
 
   std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  tt.toc();

  //  std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
  // tt.tic();
  // abs.recompute_mr(tf_new_com,tf_old, distance, w2_lb, w2_ub, rs_post);
 
  //  std::cout << "Number of new transitions: " << tf_new_com.get_no_transitions() << std::endl;
  // tt.toc();
  /* define target set */
  auto target = [&s_eta, &z, &ss](const scots::abs_type& abs_state) {
    state_type t_lb = {{0,0,0,-3*M_PI/180,-3*M_PI/180,-3*M_PI/180}};
    state_type t_ub = {{5,5,2.5,0,0,0}};
    state_type c_lb;
    state_type c_ub;
    /* center of cell associated with abs_state is stored in x */
    state_type x;
    ss.itox(abs_state,x);
    /* hyper-interval of the quantizer symbol with perturbation */
    for(int i=0; i<state_dim; i++) {
      c_lb[i] = x[i]-s_eta[i]/2.0-z[i];
      c_ub[i] = x[i]+s_eta[i]/2.0+z[i];
    }
    if( t_lb[0]<=c_lb[0] && c_ub[0]<=t_ub[0] &&
        t_lb[1]<=c_lb[1] && c_ub[1]<=t_ub[1] &&
        t_lb[2]<=c_lb[2] && c_ub[2]<=t_ub[2] &&
        t_lb[3]<=c_lb[3] && c_ub[3]<=t_ub[3] &&
        t_lb[4]<=c_lb[4] && c_ub[4]<=t_ub[4] &&
        t_lb[5]<=c_lb[5] && c_ub[5]<=t_ub[5]) {
      
        return true;
      
    }
    return false;
  };
  /* write grid point IDs with uniform grid information to file */
  write_to_file(ss,target,"target");
 
  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::solve_reachability_game(tf_new,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";

  return 1;
}
