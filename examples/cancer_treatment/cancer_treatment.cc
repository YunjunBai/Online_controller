

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

const double alphax=0.0204;
const double alphay = 0.0242;
const double betax= 0.0076;
const double betay= 0.0168;
const double k1 = 0.0;
const double k2 = 2.0;
const double k3=  0.8;
const double k4 = 0.5;
const double  m1=  0.00005;
const double z0 = 30.0;
const double t = 12.5;
const double r1 = 15.0;
const double r0 = 10.0;
const double c1 = 0.0;
const double c2 = 0.0;
const double c3 = 0.0;
const double d0 = 1.0;
const double scale =50.0;
/* state space dim */
const int state_dim=4;
/* input space dim */
const int input_dim=1;

/* sampling time */
const double tau = 12.5;

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


int main(){ 
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{0,0,0,0}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{30,30,30,30}};
  /* grid node distance diameter */
  state_type s_eta={{1,1,1,1}};
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();
  
  /* construct grid for the input alphabet */
  /* hyper-rectangle [1,2] with grid node distance 1 */
  scots::UniformGrid is(input_dim,input_type{{1}},input_type{{2}},input_type{{1}});
  is.print_info();

  

  disturbance_type w_1={{0.0005, 0.00005, 0.0005,0.0001}};
  disturbance_type w_2={{0.0003, 0.001, 0.00005,0.0001}};
  disturbance_type w2_lb={{20,20,20,20}};
  disturbance_type w2_ub={{35,25,30,30}};
  double persent=(w2_ub[0]-w2_lb[0])*(w2_ub[1]-w2_lb[1])*(w2_ub[2]-w2_lb[2])/((s_ub[0]-s_lb[0])*(s_ub[1]-s_lb[1])*(s_ub[2]-s_lb[2]));
  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);

  auto rs_post = [&dis,w2_ub,w2_lb](ds_type &y, input_type &u) -> void {
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
    double L[4][4];
    L[0][0]=5.7412;
    L[0][2]=10.7989;
    L[1][0]=0.0025;
    L[1][1]=0.3700;
    L[1][2]=1.2125;
    L[2][2]=-4*u[0];
    L[3][0]=0.5775;
    L[3][1]=0.3700;
    L[3][2]=10.7964;
    double Gx =((alphax * (k1 + ((1 - k1) * (y[2] / (y[2] + k2))))) - (betax * (k3 + (1 - k3) * (y[2] / (y[2] + k4)))));
    double Gy  =((alphay * (1 - (d0 * (y[2] / z0)))) - betay);
    double Mxy =(m1 * (1 - (y[2] / z0)));

    yy[0] = scale * (((Gx - Mxy) * y[0]) + c1 * y[0]);
    yy[1] = scale * (((Mxy * y[0]) + Gy * y[1] ) + c2 * y[1]);
    if(u[0]==0)
      yy[2] = scale * (((z0 - y[2]) / t) + c3 * y[2]);
    else
      yy[2]= scale * (((0 - y[2]) / t) + c3 * y[2]);
    yy[3] = scale * ((((Gx - Mxy) * y[0]) + c1 * y[0]) + (((Mxy * y[0]) + Gy * y[1] ) + c2 * y[1]));
    yy[4] = L[0][0]*y[4]+L[0][2]*y[6]+w[0];
    yy[5] = L[1][0]*y[4]+L[1][1]*y[5]+L[1][2]*y[6]+w[1];
    yy[6] = L[2][2]*y[6]+w[2];
    yy[7] = L[3][0]*y[4]+L[3][1]*y[5]+L[3][2]*y[6]+w[3];
   
  };
  //while(ignore==false){
    scots::runge_kutta_fixed4(rhs,y,u,dis, w2_lb,w2_ub, 2*state_dim,tau,10);
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
    //disturbance_type w = {0.05, 0.05, 0.05};
   double L[4][4];
    L[0][0]=5.7412;
    L[0][2]=10.7989;
    L[1][0]=0.0025;
    L[1][1]=0.3700;
    L[1][2]=1.2125;
    L[2][2]=-4*u[0];
    L[3][0]=0.5775;
    L[3][1]=0.3700;
    L[3][2]=10.7964;
    double Gx =((alphax * (k1 + ((1 - k1) * (y[2] / (y[2] + k2))))) - (betax * (k3 + (1 - k3) * (y[2] / (y[2] + k4)))));
    double Gy  =((alphay * (1 - (d0 * (y[2] / z0)))) - betay);
    double Mxy =(m1 * (1 - (y[2] / z0)));

    yy[0] = scale * (((Gx - Mxy) * y[0]) + c1 * y[0]);
    yy[1] = scale * (((Mxy * y[0]) + Gy * y[1] ) + c2 * y[1]);
    if(u[0]==0)
      yy[2] = scale * (((z0 - y[2]) / t) + c3 * y[2]);
    else
      yy[2]= scale * (((0 - y[2]) / t) + c3 * y[2]);
    yy[3] = scale * ((((Gx - Mxy) * y[0]) + c1 * y[0]) + (((Mxy * y[0]) + Gy * y[1] ) + c2 * y[1]));
    yy[4] = L[0][0]*y[4]+L[0][2]*y[6]+w[0];
    yy[5] = L[1][0]*y[4]+L[1][1]*y[5]+L[1][2]*y[6]+w[1];
    yy[6] = L[2][2]*y[6]+w[2];
    yy[7] =L[3][0]*y[4]+L[3][1]*y[5]+L[3][2]*y[6]+w[3];
   
  };
  //ignore = dis.get_out_of_domain();
  //if(ignore==false)    
  scots::runge_kutta_fixed4(rhs,y,u,dis,w2_lb,w2_ub, 2*state_dim,tau,10);

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
  abs.compute_gb(tf_o1d,rs_post);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();

  //if(!getrusage(RUSAGE_SELF, &usage))
 //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_o1d.get_no_transitions() << std::endl;
  dis.update_disturbance(w_2, w2_lb, w2_ub);

   std::cout << "\nComputing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post);
  
 if(!getrusage(RUSAGE_SELF, &usage))
   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  double t2= tt.toc();
  
  std::cout << "\nComputing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new,tf_o1d, w2_lb, w2_ub, rs_repost);
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  double t1=tt.toc();

  /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);

    /* function returns 1 if cell associated with x is in target set  */
    if (0 <= (x[0]-s_eta[0]/2.0) && 0<=(x[1]+s_eta[1]/2.0)  && 
        0 <= (x[2]-s_eta[2]/2.0) && 0<=(x[3]+s_eta[3]/2.0))
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
  
  return 1;
}
