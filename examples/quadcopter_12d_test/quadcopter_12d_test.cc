/*
 * quadcopter.cc
 *
 *  created: 1 2019
 *  author: Yunjun
 *
 */

/*
 * information about this example is given in
 * https://www.kth.se/polopoly_fs/1.588039.1441112632!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf
 *http://sal.aalto.fi/publications/pdf-files/eluu11_public.pdf
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
const double tau = 0.01;
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
  state_type s_eta={{2,2,2,40*M_PI/180,40*M_PI/180,40*M_PI/180,2,2,2,10*M_PI/180,10*M_PI/180,10*M_PI/180}}; 
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{-2.5,-2.5,-2.5,-20*M_PI/180,-20*M_PI/180,-20*M_PI/180,-2,-2,-2,-10*M_PI/180,-10*M_PI/180,-10*M_PI/180}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{2,2,2,60*M_PI/180,60*M_PI/180,60*M_PI/180,2,2,2,10*M_PI/180,10*M_PI/180,10*M_PI/180}}; 
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();

  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{550,550,550,550}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{675,675,675,675}};
  /* grid node distance diameter */
  input_type i_eta={{50,50,50,50}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* setup object to compute the transition function */
  scots::Abstraction<state_type,input_type,ds_type> abs(ss,is);
  
  disturbance_type w_1={.108,0.2,0,0,0,0,0.1,0,0,0,0,0};
  disturbance_type w_2={0.203, 0.1, 0.1,0,0,0,0,0,0,0,0,0};
  disturbance_type w2_lb={0,0,0,-20*M_PI/180,-20*M_PI/180,-20*M_PI/180,0,0,0,0,0,0};
  disturbance_type w2_ub={1,1,1,0,0,0,1,1,1,10*M_PI/180,10*M_PI/180,10*M_PI/180};

  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);
  

   auto l_matrix=[&is](const abs_type& input_id){
    matrix_type L;
    input_type u;
    double b=2.980*1e-6;
    double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    is.itox(input_id,u);
    L[0][3]=1.7597;
    L[0][4]= 1.3624;
    L[0][5]=1.3616;
    L[0][6]=0.9848;
    L[0][7]= 0.9962;
    L[0][8]= 0.9627;
    L[1][3]=1.7597;
    L[1][4]=1.3624;
    L[1][5]=1.3616;
    L[1][6]=0.9848;
    L[1][7]=0.9962;
    L[1][8]=0.9627;
    L[2][4]= -0.0911;
    L[2][5]=0.9698;
    L[2][6]= -0.1736;
    L[2][7]=0.8529;
    L[2][8]= 0.9698;
    L[3][4]=0.8550;
    L[3][5]=0.3438;
    L[3][10]=1.7321;
    L[3][11]=1.9696;
    L[4][5]=-3.1948e-06;
    L[4][10]= 0.9848;
    L[4][11]=-0.1736;
    L[5][4]=0.9873;
    L[5][5]=0.2977;
    L[5][9]=1;
    L[5][10]=1.5;
    L[5][11]=1.7057;
    L[6][3]=1.3715*f_t;
    L[6][4]=1.6602*f_t;
    L[6][5]=2.1287*f_t;
    L[7][3]=2.0570*f_t;
    L[7][4]=1.7438*f_t;
    L[7][5]=2.0388*f_t;
    L[8][4]= 1.8224*f_t;
    L[8][5]=1.8224*f_t;
    L[9][10]= -2.0032e-08;
    L[9][11]= -2.0032e-08;
    L[10][9]=0.1418;
    L[10][11]=0.1418;
    L[11][9]=0;
    L[11][10]=0;
   
    return L;
  };

  scots::GbEstimation<disturbance_type,matrix_type> ge(is, ss,w_1,w_2);
  ge.exp_interals(l_matrix,tau/10);

  auto rs_post = [&dis,&ge,w2_lb,w2_ub](ds_type &y, input_type &u) -> void {
  auto rhs_1 =[&dis, &ge](ds_type &yy, const ds_type &y, input_type &u) -> void {
    /* find the distrubance for the given state */
    state_type x;
    state_type r;
    state_type r_es;
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    //std::cout<<"test test x:"<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<" "<<x[4]<<" "<<x[5]<<" "<<x[6]<<" "<<x[7]<<" "<<x[8]<<" "<<x[9]<<" "<<x[10]<<" "<<x[11]<<std::endl;
    //std::cout<<"test test r:"<<r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<r[3]<<" "<<r[4]<<" "<<r[5]<<" "<<r[6]<<" "<<r[7]<<" "<<r[8]<<" "<<r[9]<<" "<<r[10]<<" "<<r[11]<<std::endl;
      r_es=ge.gb_estimate(r,u);
    disturbance_type w = dis.get_disturbance(x,r_es); 
   
    double m=0.468;
    double g=9.81;
    double b=2.980*1e-6;
    double l=0.225;
    double d=1.140*1e-7;
    double I_x=4.856*1e-3;
    double I_y=4.856*1e-3;
    double I_z=8.801*1e-3;
    double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    double tau_x=b*l*(u[2]*u[2]-u[0]*u[0]);
    double tau_y=b*l*(u[3]*u[3]-u[1]*u[1]);
    double tau_z=d*(u[3]*u[3]+u[1]*u[1]-u[2]*u[2]-u[0]*u[0]);
     double L[12][12];
    L[0][3]=1.7597;
    L[0][4]= 1.3624;
    L[0][5]=1.3616;
    L[0][6]=0.9848;
    L[0][7]= 0.9962;
    L[0][8]= 0.9627;
    L[1][3]=1.7597;
    L[1][4]=1.3624;
    L[1][5]=1.3616;
    L[1][6]=0.9848;
    L[1][7]=0.9962;
    L[1][8]=0.9627;
    L[2][4]= -0.0911;
    L[2][5]=0.9698;
    L[2][6]= -0.1736;
    L[2][7]=0.8529;
    L[2][8]= 0.9698;
    L[3][4]=0.8550;
    L[3][5]=0.3438;
    L[3][10]=1.7321;
    L[3][11]=1.9696;
    L[4][5]=-3.1948e-06;
    L[4][10]=  0.9848;
    L[4][11]=-0.1736;
    L[5][4]=0.9873;
    L[5][5]=0.2977;
    L[5][9]=1;
    L[5][10]=1.5;
    L[5][11]=1.7057;
    L[6][3]=1.3715*f_t;
    L[6][4]=1.6602*f_t;
    L[6][5]=2.1287*f_t;
    L[7][3]=2.0570*f_t;
    L[7][4]=1.7438*f_t;
    L[7][5]=2.0388*f_t;
    L[8][4]= 1.8224*f_t;
    L[8][5]=1.8224*f_t;
    L[9][10]= -2.0032e-08;
    L[9][11]= -2.0032e-08;
    L[10][9]=0.1418;
    L[10][11]=0.1418;
    L[11][9]=0;
    L[11][10]=0;

    yy[0] = y[8]*(std::sin(y[5])*std::sin(y[3])+std::cos(y[5])*std::cos(y[3])*std::sin(y[4])) - y[7]*(std::cos(y[5])*std::sin(y[3])-std::cos(y[3])*std::sin(y[5])*std::sin(y[4])) +y[6]*std::cos(y[3])*std::cos(y[4]);
    yy[1] =-y[8]*(std::cos(y[3])*std::sin(y[5])-std::cos(y[5])*std::sin(y[3])*std::sin(y[4])) + y[7]*(std::cos(y[5])*std::cos(y[3])+std::sin(y[5])*std::sin(y[3])*std::sin(y[4])) +y[6]*std::cos(y[4])*std::sin(y[3]);
    yy[2] =y[8]*std::cos(y[5])*std::cos(y[4])-y[6]*std::sin(y[4])+y[7]*std::cos(y[4])*std::sin(y[5]);
    yy[3] =y[10]*(std::sin(y[5])/std::cos(y[4]))+y[11]*(std::cos(y[5])/std::cos(y[4]));
    yy[4] =y[10]*std::cos(y[5])-y[11]*std::sin(y[5]);
    yy[5] =y[9]+y[10]*std::sin(y[5])*std::tan(y[4])+y[11]*std::cos(y[5])*std::tan(y[4]);
    yy[6] =1/m *(std::sin(y[5])*std::sin(y[3])+std::cos(y[5])*std::cos(y[3])*std::sin(y[4])) * f_t;
    yy[7] =1/m *(std::cos(y[3])*std::sin(y[5])-std::cos(y[5])*std::sin(y[3])*std::sin(y[4])) *f_t;
    yy[8] =1/m *std::cos(y[5])*std::cos(y[4])*f_t+g;
    yy[9] =(I_y-I_z)/I_x*y[10]*y[11] + 1/I_x*tau_x;
    yy[10]=(I_z-I_x)/I_y*y[9]*y[11] + 1/I_y*tau_y;
    yy[11]=(I_x-I_y)/I_z*y[9]*y[10] +1/I_z *tau_z;

        /* to account for input disturbances */
    yy[12] = L[0][3]*y[15]+L[0][4]*y[16]+L[0][5]*y[17]+L[0][6]*y[18]+L[0][7]*y[19]+L[0][8]*y[20]+w[0];
    yy[13] = L[1][3]*y[15]+L[1][4]*y[16]+L[1][5]*y[17]+L[1][6]*y[18]+L[1][7]*y[19]+L[1][8]*y[20]+w[1];
    yy[14] = L[2][4]*y[16]+L[2][5]*y[17]+L[2][6]*y[18]+L[2][7]*y[19]+L[2][8]*y[20]+w[2];
    yy[15] = L[3][4]*y[16]+L[3][5]*y[17]+L[3][10]*y[22]+L[3][11]*y[23]+w[3];
    yy[16] = L[4][5]*y[17]+L[4][10]*y[22]+L[4][11]*y[23]+w[4];
    yy[17] = L[5][4]*y[16]+L[5][5]*y[17]+L[5][9]*y[21]+L[5][10]*y[22]+L[5][11]*y[23]+w[5];
    yy[18] = L[6][3]*y[15]+L[6][4]*y[16]+L[6][5]*y[17]+w[6];
    yy[19] = L[7][3]*y[15]+L[7][4]*y[16]+L[7][5]*y[17]+w[7];
    yy[20] = L[8][4]*y[16]+L[8][5]*y[17]+w[8];
    yy[21] = L[9][10]*y[22]+L[9][11]*y[23]+w[9];
    yy[22] = L[10][9]*y[21]+L[9][11]*y[23]+w[10];
    yy[23] = L[11][9]*y[21]+L[11][10]*y[22]+w[11];
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
    for (int i=0; i<state_dim; i++){
      x[i] = y[i];
      r[i] = y[i+state_dim];
    }
    
    r_es=ge.gb_estimate(r,u);
    disturbance_type w = dis.get_disturbance(x,r_es); 
     double m=0.468;
    double g=9.81;
    double b=2.980*1e-6;
    double l=0.225;
    double d=1.140*1e-7;
    double I_x=4.856*1e-3;
    double I_y=4.856*1e-3;
    double I_z=8.801*1e-3;
    double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    double tau_x=b*l*(u[2]*u[2]-u[0]*u[0]);
    double tau_y=b*l*(u[3]*u[3]-u[1]*u[1]);
    double tau_z=d*(u[3]*u[3]+u[1]*u[1]-u[2]*u[2]-u[0]*u[0]);
     double L[12][12];
    L[0][3]=1.7597;
    L[0][4]= 1.3624;
    L[0][5]=1.3616;
    L[0][6]=0.9848;
    L[0][7]= 0.9962;
    L[0][8]= 0.9627;
    L[1][3]=1.7597;
    L[1][4]=1.3624;
    L[1][5]=1.3616;
    L[1][6]=0.9848;
    L[1][7]=0.9962;
    L[1][8]=0.9627;
    L[2][4]= -0.0911;
    L[2][5]=0.9698;
    L[2][6]= -0.1736;
    L[2][7]=0.8529;
    L[2][8]= 0.9698;
    L[3][4]=0.8550;
    L[3][5]=0.3438;
    L[3][10]=1.7321;
    L[3][11]=1.9696;
    L[4][5]=-3.1948e-06;
    L[4][10]=  0.9848;
    L[4][11]=-0.1736;
    L[5][4]=0.9873;
    L[5][5]=0.2977;
    L[5][9]=1;
    L[5][10]=1.5;
    L[5][11]=1.7057;
    L[6][3]=1.3715*f_t;
    L[6][4]=1.6602*f_t;
    L[6][5]=2.1287*f_t;
    L[7][3]=2.0570*f_t;
    L[7][4]=1.7438*f_t;
    L[7][5]=2.0388*f_t;
    L[8][4]= 1.8224*f_t;
    L[8][5]=1.8224*f_t;
    L[9][10]= -2.0032e-08;
    L[9][11]= -2.0032e-08;
    L[10][9]=0.1418;
    L[10][11]=0.1418;
    L[11][9]=0;
    L[11][10]=0;
    
    yy[0] = y[8]*(std::sin(y[5])*std::sin(y[3])+std::cos(y[5])*std::cos(y[3])*std::sin(y[4])) - y[7]*(std::cos(y[5])*std::sin(y[3])-std::cos(y[3])*std::sin(y[5])*std::sin(y[4])) +y[6]*std::cos(y[3])*std::cos(y[4]);
    yy[1] =-y[8]*(std::cos(y[3])*std::sin(y[5])-std::cos(y[5])*std::sin(y[3])*std::sin(y[4])) + y[7]*(std::cos(y[5])*std::cos(y[3])+std::sin(y[5])*std::sin(y[3])*std::sin(y[4])) +y[6]*std::cos(y[4])*std::sin(y[3]);
    yy[2] =y[8]*std::cos(y[5])*std::cos(y[4])-y[6]*std::sin(y[4])+y[7]*std::cos(y[4])*std::sin(y[5]);
    yy[3] =y[10]*(std::sin(y[5])/std::cos(y[4]))+y[11]*(std::cos(y[5])/std::cos(y[4]));
    yy[4] =y[10]*std::cos(y[5])-y[11]*std::sin(y[5]);
    yy[5] =y[9]+y[10]*std::sin(y[5])*std::tan(y[4])+y[11]*std::cos(y[5])*std::tan(y[4]);
    yy[6] =-1/m *(std::sin(y[5])*std::sin(y[3])+std::cos(y[5])*std::cos(y[3])*std::sin(y[4])) * f_t;
    yy[7] =-1/m *(std::cos(y[3])*std::sin(y[5])-std::cos(y[5])*std::sin(y[3])*std::sin(y[4])) *f_t;
    yy[8] =-1/m *std::cos(y[5])*std::cos(y[4])*f_t+g;
    yy[9] =(I_y-I_z)/I_x*y[10]*y[11] + 1/I_x*tau_x;
    yy[10]=(I_z-I_x)/I_y*y[9]*y[11] + 1/I_y*tau_y;
    yy[11]=(I_x-I_y)/I_z*y[9]*y[10] +1/I_z *tau_z;

        /* to account for input disturbances */
    yy[12] = L[0][3]*y[15]+L[0][4]*y[16]+L[0][5]*y[17]+L[0][6]*y[18]+L[0][7]*y[19]+L[0][8]*y[20]+w[0];
    yy[13] = L[1][3]*y[15]+L[1][4]*y[16]+L[1][5]*y[17]+L[1][6]*y[18]+L[1][7]*y[19]+L[1][8]*y[20]+w[1];
    yy[14] = L[2][4]*y[16]+L[2][5]*y[17]+L[2][6]*y[18]+L[2][7]*y[19]+L[2][8]*y[20]+w[2];
    yy[15] = L[3][4]*y[16]+L[3][5]*y[17]+L[3][10]*y[22]+L[3][11]*y[23]+w[3];
    yy[16] = L[4][5]*y[17]+L[4][10]*y[22]+L[4][11]*y[23]+w[4];
    yy[17] = L[5][4]*y[16]+L[5][5]*y[17]+L[5][9]*y[21]+L[5][10]*y[22]+L[5][11]*y[23]+w[5];
    yy[18] = L[6][3]*y[15]+L[6][4]*y[16]+L[6][5]*y[17]+w[6];
    yy[19] = L[7][3]*y[15]+L[7][4]*y[16]+L[7][5]*y[17]+w[7];
    yy[20] = L[8][4]*y[16]+L[8][5]*y[17]+w[8];
    yy[21] = L[9][10]*y[22]+L[9][11]*y[23]+w[9];
    yy[22] = L[10][9]*y[21]+L[9][11]*y[23]+w[10];
    yy[23] = L[11][9]*y[21]+L[11][10]*y[22]+w[11];
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
  std::queue<abs_type> online_queue; 

  std::cout << "Computing the transition function: " << std::endl;
  tt.tic();
  abs.compute_gb(tf_old, rs_post);
  tt.toc();
   if(!getrusage(RUSAGE_SELF, &usage))
     std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_old.get_no_transitions() << std::endl;
   std::cout << "Number of transitions: " << tf_old.get_no_transitions() << std::endl;

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
  abs.recompute_gb(tf_new,online_queue,tf_old, w2_lb, w2_ub, rs_repost);
 
   std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  tt.toc();

  //  std::cout << "Computing the new transition function locally (after distrubance changes): " << std::endl;
  // tt.tic();
  // abs.recompute_mr(tf_new_com,tf_old, distance, w2_lb, w2_ub, rs_post);
 
  //  std::cout << "Number of new transitions: " << tf_new_com.get_no_transitions() << std::endl;
  // tt.toc();
  /* define target set */
  auto target = [&s_eta, &ss](const scots::abs_type& abs_state) {
    state_type t_lb = {{0.6,0.6,0.6,-20*M_PI/180,-20*M_PI/180,-20*M_PI/180,0,0,0,0,0,0}};
    state_type t_ub = {{1,1,1,-20*M_PI/180,-20*M_PI/180,-20*M_PI/180,0,0,0,0,0,0}};
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
