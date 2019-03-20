/*
 * quadcopter.cc
 *
 *  created: 1 2019
 *  author: Yunjun
 *
 */

/*
 * information about this example is given in
 "Control and Estimation of a Quadcopter Dynamical Model"
Sevkuthan KURAK1

FENS, Inernational University of Sarajevo 
 * https://www.kth.se/polopoly_fs/1.588039.1441112632!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf
 *http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=8847641&fileOId=8859343
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
const int state_dim=12;
/* input space dim */
const int input_dim=4;
/* sampling time */
const double tau = 0.001;
using abs_type = scots::abs_type;
/* data types of the state space elements and input 
 * space elements used in uniform grid and ode solver */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;
using disturbance_type = std::array<double, state_dim>;
using ds_type = std::array<double, 2*state_dim>;
using matrix_type = std::array<std::array<double,state_dim>,state_dim>;

int main() {
  /* to measure time */
  TicToc tt;

  /* construct grid for the state space */
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  state_type s_eta={{1.5,1.5,1.5,20*M_PI/180,20*M_PI/180,20*M_PI/180,0.3,0.3,0.3,0.05,0.05,0.05}}; 
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{-2.5,-2.5,-2.5,-30*M_PI/180,-30*M_PI/180,-30*M_PI/180,0,0,0,-0.05,-0.05,-0.05}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{2,2,2,30*M_PI/180,30*M_PI/180,30*M_PI/180,0.9,0.9,0.9,0.05,0.05,0.05}}; 
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();

  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{100,100,100,100}};
  /* upper bounds of the hyper rectangle */
  //input_type i_ub={{2500,2500,2500,2500}};
  input_type i_ub={{300,300,300,300}};
  /* grid node distance diameter */
  input_type i_eta={{50,50,50,50}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* setup object to compute the transition function */
  scots::Abstraction<state_dim,input_dim> abs(ss,is);
   disturbance_type w_1={0, 0, 0,0,0,0,0,0,0,0.123,0.123,0};
  disturbance_type w_2={0, 0, 0,0,0,0,0,0,0,.503,.503,0};
  disturbance_type w2_lb={-2.5,-2.5,-2.5,-30*M_PI/180,-30*M_PI/180,-30*M_PI/180,0,0,0,-0.05,-0.05,-0.05};
  disturbance_type w2_ub={2,2,2,30*M_PI/180,30*M_PI/180,30*M_PI/180,0.4,0.9,0.9,0.05,0.05,0.05};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);

   auto l_matrix=[&is](const abs_type& input_id){
    matrix_type L;
    input_type u;
    double b=54.2 * 1e-6;
    double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    is.itox(input_id,u);
    L[0][6]=1;
    L[1][7]=1;
    L[2][8]= 1;
    L[3][9]=1;
    L[4][10]=1;
    L[5][11]=1;
    L[6][3]=1.3715*f_t;
    L[6][4]=1.6602*f_t;
    L[6][5]=2.1287*f_t;
    L[7][3]=2.0570*f_t;
    L[7][4]=1.7438*f_t;
    L[7][5]=2.0388*f_t;
    L[8][4]= 1.8224*f_t;
    L[8][5]=1.8224*f_t;
    L[9][10]= 2.0032e-08;
    L[9][11]= 2.0032e-08;
    L[10][9]=0.1418;
    L[10][11]=0.1418;
   
    return L;
  };

  scots::GbEstimation<disturbance_type,matrix_type> ge(is, ss,w_1,w_2);
  ge.exp_interals(l_matrix,tau/10);
  int test_counter=10;
  auto rs_post = [&dis,&ge,w2_lb,w2_ub,&test_counter](ds_type &y, input_type &u) -> void {
  auto rhs_1 =[&dis, &ge,&test_counter](ds_type &yy, const ds_type &y, input_type &u) -> void {
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
     w = dis.get_disturbance(x,r_es);

    double m=0.506;
    double g=9.81;
    double b=54.2 * 1e-6;
    double l=0.24;
    double d=1.1*1e-6;
    double I_x=8.1*1e-5;
    double I_y=8.1*1e-5;
    double I_z=6.2*1e-5;
    double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    double tau_x=b*l*(u[2]*u[2]-u[0]*u[0]);
    double tau_y=b*l*(u[3]*u[3]-u[1]*u[1]);
    double tau_z=d*(u[3]*u[3]+u[1]*u[1]-u[2]*u[2]-u[0]*u[0]);
     double L[12][12];
    
    L[0][6]=1;
    L[1][7]=1;
    L[2][8]= 1;
    L[3][9]=1;
    L[4][10]=1;
    L[5][11]=1;
    L[6][3]=1.3715*f_t;
    L[6][4]=1.6602*f_t;
    L[6][5]=2.1287*f_t;
    L[7][3]=2.0570*f_t;
    L[7][4]=1.7438*f_t;
    L[7][5]=2.0388*f_t;
    L[8][4]= 1.8224*f_t;
    L[8][5]=1.8224*f_t;
    L[9][10]= 2.0032e-08;
    L[9][11]= 2.0032e-08;
    L[10][9]=0.1418;
    L[10][11]=0.1418;

    yy[0] = y[6];
    yy[1] = y[7];
    yy[2] = y[8];
    yy[3] =std::cos(y[5])/std::sin(y[4]) * y[10] +std::cos(y[5])/std::cos(y[4]) * y[11];
    yy[4] =std::cos(y[5])* y[10]+std::sin(y[5])* y[11];
    yy[5] =y[9]+y[10]*std::sin(y[5])*std::tan(y[4])+y[11]*std::cos(y[5])*std::tan(y[4]);
    yy[6] =-1/m *(std::sin(y[5])*std::sin(y[3])+std::cos(y[5])*std::cos(y[3])*std::sin(y[4])) * f_t;
    yy[7] =-1/m *(std::cos(y[3])*std::sin(y[5])-std::cos(y[5])*std::sin(y[3])*std::sin(y[4])) *f_t;
    yy[8] =-1/m *std::cos(y[5])*std::cos(y[4])*f_t+g;
    yy[9] =((I_y-I_z)/I_x)*y[10]*y[11] + 1/I_x*tau_x;
    yy[10]=((I_z-I_x)/I_y)*y[9]*y[11] + 1/I_y*tau_y;
    yy[11]=((I_x-I_y)/I_z)*y[9]*y[10] +1/I_z *tau_z;

        /* to account for input disturbances */
    yy[12] = L[0][6]*y[18]+w[0];
    yy[13] = L[1][7]*y[19]+w[1];
    yy[14] = L[2][8]*y[20]+w[2];
    yy[15] = L[3][9]*y[21]+w[3];
    yy[16] = L[4][10]*y[22]+w[4];
    yy[17] = L[5][11]*y[23]+w[5];
    yy[18] = L[6][3]*y[15]+L[6][4]*y[16]+L[6][5]*y[17]+w[6];
    yy[19] = L[7][3]*y[15]+L[7][4]*y[16]+L[7][5]*y[17]+w[7];
    yy[20] = L[8][4]*y[16]+L[8][5]*y[17]+w[8];
    yy[21] = L[9][10]*y[22]+L[9][11]*y[23]+w[9];
    yy[22] = L[10][9]*y[21]+L[9][11]*y[23]+w[10];
    yy[23] = w[11];
  };
  scots::runge_kutta_fixed4(rhs_1,y,u,dis, w2_lb,w2_ub,2*state_dim,tau,10);
};

auto rs_repost = [&dis, &ge,w2_lb,w2_ub](ds_type &y, input_type &u, bool &neigbour) -> void {
  dis.set_intersection_check();
  //dis.set_out_of_domain();
  auto rhs =[&dis,&ge,w2_lb,w2_ub](ds_type &yy, const ds_type &y, input_type &u) -> void {
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
     w = dis.get_disturbance(x,r_es);

    double m=0.506;
    double g=9.81;
    double b=54.2 * 1e-6;
    double l=0.24;
    double d=1.1*1e-6;
    double I_x=8.1*1e-5;
    double I_y=8.1*1e-5;
    double I_z=6.2*1e-5;
    double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    double tau_x=b*l*(u[2]*u[2]-u[0]*u[0]);
    double tau_y=b*l*(u[3]*u[3]-u[1]*u[1]);
    double tau_z=d*(u[3]*u[3]+u[1]*u[1]-u[2]*u[2]-u[0]*u[0]);
     double L[12][12];
    
    L[0][6]=1;
    L[1][7]=1;
    L[2][8]= 1;
    L[3][9]=1;
    L[4][10]=1;
    L[5][11]=1;
    L[6][3]=1.3715*f_t;
    L[6][4]=1.6602*f_t;
    L[6][5]=2.1287*f_t;
    L[7][3]=2.0570*f_t;
    L[7][4]=1.7438*f_t;
    L[7][5]=2.0388*f_t;
    L[8][4]= 1.8224*f_t;
    L[8][5]=1.8224*f_t;
    L[9][10]= 2.0032e-08;
    L[9][11]= 2.0032e-08;
    L[10][9]=0.1418;
    L[10][11]=0.1418;

    yy[0] = y[6];
    yy[1] = y[7];
    yy[2] = y[8];
    yy[3] =std::cos(y[5])/std::sin(y[4]) * y[10] +std::cos(y[5])/std::cos(y[4]) * y[11];
    yy[4] =std::cos(y[5])* y[10]+std::sin(y[5])* y[11];
    yy[5] =y[9]+y[10]*std::sin(y[5])*std::tan(y[4])+y[11]*std::cos(y[5])*std::tan(y[4]);
    yy[6] =-1/m *(std::sin(y[5])*std::sin(y[3])+std::cos(y[5])*std::cos(y[3])*std::sin(y[4])) * f_t;
    yy[7] =-1/m *(std::cos(y[3])*std::sin(y[5])-std::cos(y[5])*std::sin(y[3])*std::sin(y[4])) *f_t;
    yy[8] =-1/m *std::cos(y[5])*std::cos(y[4])*f_t+g;
    yy[9] =((I_y-I_z)/I_x)*y[10]*y[11] + 1/I_x*tau_x;
    yy[10]=((I_z-I_x)/I_y)*y[9]*y[11] + 1/I_y*tau_y;
    yy[11]=((I_x-I_y)/I_z)*y[9]*y[10] +1/I_z *tau_z;

        /* to account for input disturbances */
    yy[12] = L[0][6]*y[18]+w[0];
    yy[13] = L[1][7]*y[19]+w[1];
    yy[14] = L[2][8]*y[20]+w[2];
    yy[15] = L[3][9]*y[21]+w[3];
    yy[16] = L[4][10]*y[22]+w[4];
    yy[17] = L[5][11]*y[23]+w[5];
    yy[18] = L[6][3]*y[15]+L[6][4]*y[16]+L[6][5]*y[17]+w[6];
    yy[19] = L[7][3]*y[15]+L[7][4]*y[16]+L[7][5]*y[17]+w[7];
    yy[20] = L[8][4]*y[16]+L[8][5]*y[17]+w[8];
    yy[21] = L[9][10]*y[22]+L[9][11]*y[23]+w[9];
    yy[22] = L[10][9]*y[21]+L[9][11]*y[23]+w[10];
    yy[23] = w[11];
  };
  //ignore = dis.get_out_of_domain();
  //if(ignore==false)    
  scots::runge_kutta_fixed4(rhs,y,u,dis, w2_lb,w2_ub,2*state_dim,tau,10);

  if(dis.get_intersection_check()==true){
    neigbour=true;
  }
};
//disturbance_type max_distance={2*tau,2*tau,2*tau,4*tau,4*tau,6*tau,2*tau,2*tau,2*tau,2*tau,2*tau,2*tau};
//disturbance_type max_w2_lb={1-2*tau,1-2*tau,1-2*tau,-30*M_PI/180,-30*M_PI/180,-30*M_PI/180,0,0,0,-0.05,-0.05,-0.05};
//disturbance_type max_w2_ub={1.5+2*tau,1.5+2*tau,1.5+2*tau,10*M_PI/180,10*M_PI/180,10*M_PI/180+6*tau,0.4+2*tau,0.4+2*tau,0.4+2*tau,0.05+2*tau,0.05+2*tau,0.05+2*tau};

 /* transition function of symbolic model */
  scots::TransitionFunction tf_old, tf_new, tf_standard,tf_new_com;
   std::queue<abs_type> online_queue; 

  std::cout << "Computing the transition function: " << std::endl;
  tt.tic();
  abs.compute_gb(tf_old, rs_post);
  tt.toc();
   if(!getrusage(RUSAGE_SELF, &usage))
     std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_old.get_no_transitions() << std::endl;
   std::cout << "Number of transitions: " << tf_old.get_no_transitions() << std::endl;

 
      dis.update_disturbance(w_2, w2_lb, w2_ub);
      std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
      tt.tic();
      abs.recompute_gb(tf_new,online_queue,tf_old, w2_lb, w2_ub, rs_repost);
 
   std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  tt.toc();

 std::cout << "Computing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post);
  
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_standard.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  tt.toc();

  /* define target set */
  auto target = [&s_eta, &ss](const scots::abs_type& abs_state) {
    state_type t_lb = {{0,0,-2.5,-30*M_PI/180,-30*M_PI/180,-30*M_PI/180,0,0,0,-0.05,-0.05,-0.05}};
    state_type t_ub = {{2,2,2,30*M_PI/180,30*M_PI/180,30*M_PI/180,0.9,0.9,0.9,0.05,0.05,0.05}};
   
    state_type c_lb;
    state_type c_ub;
    /* center of cell associated with abs_state is stored in x */
    state_type x;
    ss.itox(abs_state,x);
    /* hyper-interval of the quantizer symbol with perturbation */
    for(int i=0; i<state_dim; i++) {
      c_lb[i] = x[i]-s_eta[i]/2.0;
      c_ub[i] = x[i]+s_eta[i]/2.0;
    }
    if( t_lb[0]<=c_lb[0] && c_ub[0]<=t_ub[0] &&
        t_lb[1]<=c_lb[1] && c_ub[1]<=t_ub[1] ) {
      
        return true;
      
    }
    return false;
  };
  /* write grid point IDs with uniform grid information to file */
  write_to_file(ss,target,"target");
 
  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::solve_reachability_game(tf_old,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
    std::cout << "Done. \n";
   return 1;
}


