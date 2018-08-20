/*
 * dcdc.cc
 *
 *  created: Oct 2015
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in the readme file
 */

#include <iostream>
#include <array>
#include <cmath>

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
const int input_dim=1;
/* sampling time */
const double tau = 0.00000001;

/*
 * data types for the elements of the state space 
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using disturbance_type = std::array<double, state_dim>;
using input_type = std::array<double,input_dim>;
using ds_type = std::array<double, 2*state_dim>;
/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* parameters for system dynamics */
const double xc=70;
const double xl=3;
const double rc=0.005;
const double rl=0.05;
const double ro=1;
const double vs=1;
/* parameters for radius calculation */
const double mu=std::sqrt(2);
/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */


int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{2.0/4e3,2.0/4e3}};
 /* lower bounds of the hyper-rectangle */
  state_type lb={{1.15,5.45}};
  /* upper bounds of the hyper-rectangle */
  state_type ub={{1.55,5.85}};
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  std::cout << "Uniform grid details:\n";
  ss.print_info();

  /* construct grid for the input alphabet */
  /* hyper-rectangle [1,2] with grid node distance 1 */
  scots::UniformGrid is(input_dim,input_type{{1}},input_type{{2}},input_type{{1}});
  is.print_info();

  disturbance_type w_1={{1/3, 0}};
  disturbance_type w_2={{0.8/3, 0}};
  disturbance_type w2_lb={{1.15,5.45}};
  disturbance_type w2_ub={{1.30,5.55}};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);

  auto rs_post = [&dis](ds_type &y, input_type &u) -> void {
  auto rhs =[&dis](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    disturbance_type w = dis.get_disturbance(x,r);

      if(u[0]==1) {
      yy[0]=-rl/xl*y[0]+w[0];
      yy[1]=-1/(xc*(ro+rc))*y[1]+w[1];
      yy[2]=-rl/xl*y[2];
      yy[3]=-1/(xc*(ro+rc))*y[3];
    } else {
      yy[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*y[0]-(1/xl)*ro/(5*(ro+rc))*y[1]+w[0];
      yy[1]=(1/xc)*5*ro/(ro+rc)*y[0]-(1/xc)*(1/(ro+rc))*y[1];
      yy[2]=-(1/xl)*(rl+ro*rc/(ro+rc))*y[2]+(1/xl)*ro/(5*(ro+rc))*y[3];
      yy[3]=5*(1/xc)*ro/(ro+rc)*y[2]-(1/xc)*(1/(ro+rc))*y[3];
    }
  };
  scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,tau,10);
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
    if(u[0]==1) {
      yy[0]=-rl/xl*y[0]+w[0];
      yy[1]=-1/(xc*(ro+rc))*y[1]+w[1];
      yy[2]=-rl/xl*y[2];
      yy[3]=-1/(xc*(ro+rc))*y[3];
    } else {
      yy[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*y[0]-(1/xl)*ro/(5*(ro+rc))*y[1]+w[0];
      yy[1]=(1/xc)*5*ro/(ro+rc)*y[0]-(1/xc)*(1/(ro+rc))*y[1];
      yy[2]=-(1/xl)*(rl+ro*rc/(ro+rc))*y[2]+(1/xl)*ro/(5*(ro+rc))*y[3];
      yy[3]=5*(1/xc)*ro/(ro+rc)*y[2]-(1/xc)*(1/(ro+rc))*y[3];
   }
  }; 
  scots::runge_kutta_fixed4(rhs,y,u,2*state_dim,tau,10);
  if(dis.get_intersection_check()){
    neigbour=true;
  }  
};


  /* compute transition function of symbolic model */
  std::cout << "Computing the transition function:\n";
  /* transition function of symbolic model */
  scots::TransitionFunction tf_o1d,tf_new,tf_standard,tf_new_com;
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
  abs.verbose_off();

   tt.tic();
  abs.compute_gb(tf_o1d,rs_post);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();

  //if(!getrusage(RUSAGE_SELF, &usage))
 //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
 // std::cout << "Number of transitions: " << tf_o1d.get_no_transitions() << std::endl;
  dis.update_disturbance(w_2, w2_lb, w2_ub);
  state_type max_dynamic = {{ub[0],ub[1]}};
  state_type distance = dis.get_maxdistance(max_dynamic,tau);

   std::cout << "Computing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post);
  
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  tt.toc();

  std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new,tf_o1d,tf_standard, distance, w2_lb, w2_ub, rs_repost);
 
   std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  tt.toc();

   std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_mr(tf_new_com,tf_o1d, distance, w2_lb, w2_ub, rs_post);
 
   std::cout << "Number of new transitions: " << tf_new_com.get_no_transitions() << std::endl;
  tt.toc();

  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub[0] && 
        lb[1] <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= ub[1])
      return true;
    return false;
  };

/********************************************************************/






/*******************************************************************/
  /* compute winning domain (contains also valid inputs) */
  std::cout << "\nSynthesis: \n";
  tt.tic();
  scots::WinningDomain win = scots::solve_invariance_game(tf_new,safeset);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << "\n";

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";

  return 1;
}

