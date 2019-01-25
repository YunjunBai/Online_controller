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
const int state_dim=4;
/* input space dim */
const int input_dim=2;

/* sampling time */
const double tau = 0.003;
const double m_1=1;
const double m_2=1;
const double l_1=0.5;
const double l_2=0.5;
const double g=9.8;
const double A=l_1*std::pow(l_2,2) *m_2;
const double B=g*l_1*l_2*(m_1+m_2);
const double C=std::pow(l_1,2)*l_2*m_2;
const double D=g*l_1*l_2*m_2;
const double E=std::pow(l_1,2)*l_2*m_1+std::pow(l_1,2)*l_2*m_2;
const double G=g*l_1*std::pow(l_2,2)*m_2*(m_2+m_1);
const double Mm=g*std::pow(l_1,2)*l_2*m_2*(m_2+m_1);
const double Q=g*l_1*std::pow(l_2,2)*std::pow(l_2,2);
const double R=std::pow(l_1,2)*m_2+std::pow(l_1,2)*m_1+std::pow(l_2,2)*m_2;
const double Uu=l_1*l_2*m_2;
const double S=std::pow(l_1,2)*std::pow(l_2,2)*std::pow(m_2,2);
const double T=m_1*m_2*std::pow(l_1,2)*std::pow(l_2,2);

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
  state_type s_lb={{-M_PI/2,-M_PI/4,-M_PI/2,-M_PI/4}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{M_PI/2,M_PI/4,M_PI/2,M_PI/4}};
  /* grid node distance diameter */
  state_type s_eta={{M_PI/20, M_PI/40,M_PI/20,M_PI/40 }};
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Uniform grid details:" << std::endl;
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{-20,-20}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{ 20, 20}};
  /* grid node distance diameter */
  input_type i_eta={{2,2}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();


   disturbance_type w_1={{0, 0.4606, 0, 1.0742 }};
  // //disturbance_type w_1={{0, 3.4593, 0, 8.0742 }};
   disturbance_type w_2={{0, 3.4593, 0, 8.0742 }};
  // disturbance_type w2_lb={{-M_PI/2,-M_PI/4,-M_PI/2,-M_PI/4}};
  // disturbance_type w2_ub={{M_PI/2,0,M_PI/4,0}};
  //disturbance_type w_3={{0, 3.4593, 0, 8.0742 }};
  disturbance_type w3_lb={{-1.571, 0.471, 1.257, -0.157}};
  disturbance_type w3_ub={{1.885,  2.199, 4.712, 1.571}};
  // double persent=(w2_ub[0]-w2_lb[0])*(w2_ub[1]-w2_lb[1])*(w2_ub[2]-w2_lb[2])/((s_ub[0]-s_lb[0])*(s_ub[1]-s_lb[1])*(s_ub[2]-s_lb[2]));
  scots::Disturbance<disturbance_type, state_type> dis(w_1, ss);
  dis.disturbance_read();
  auto l_matrix=[&is](const abs_type& input_id){
    matrix_type l;
    input_type u;
    is.itox(input_id,u);
    l[0][1]=1;
    l[2][3]=1;
    l[1][0]=6.9296;
    l[1][1]=0.2210;
    l[1][2]=19.7559+2*(u[0]+2*u[1])-1.0005e-07*(u[0]-3*u[1])-2.3020e-07*(u[0]-2*u[1]);
    l[1][3]=0.0628;
    l[3][0]=31.8054;
    l[3][1]=0.4638;
    l[3][2]=9.8780 + 0.1113*(20*u[0]-20*u[1])+0.0369*20*u[1];
    l[3][3]=0.1013;
    return l;
  };

  scots::GbEstimation<disturbance_type,matrix_type> ge(is, ss,w_1, w_2);
  ge.exp_interals(l_matrix,tau/10);

  auto rs_post = [&dis,&ge,w3_ub,w3_lb](ds_type &y, input_type &u) -> void {
   // dis.set_out_of_domain();
  auto rhs =[&dis,&ge](ds_type &yy, const ds_type &y, input_type &u) -> void {
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
    double L[4][4];
    L[0][1]=1;
    L[2][3]=1;
    L[1][0]=6.9296;
    L[1][1]=0.2210;
    L[1][2]=19.7559+2*(u[0]+2*u[1])-1.0005e-07*(u[0]-3*u[1])-2.3020e-07*(u[0]-2*u[1]);
    L[1][3]=0.0628;
    L[3][0]=31.8054;
    L[3][1]=0.4638;
    L[3][2]=9.8780 + 0.1113*(20*u[0]-20*u[1])+0.0369*20*u[1];
    L[3][3]=0.1013;
    double F=(2*l_1*std::pow(l_2,3)*std::pow(m_2,2) + std::pow(l_1,3) *l_2*std::pow(m_2,2) + std::pow(l_1,3) *l_2*m_2*m_1)* std::pow(y[1],2)+2*l_1*std::pow(l_2,3)*std::pow(m_2,2)*y[1]*y[3];
    double H=std::pow(l_1,2)*std::pow(l_2,2)*std::pow(m_2,2)*(3*std::pow(y[2],2)+2*y[1]*y[3]);
    double P=l_1*l_2*m_2*u[0]-2*l_1*l_2*m_2*u[1];
  
    yy[0]=y[1];
    yy[1]=(l_2*(u[0]-u[1])+(D*std::sin(y[0]+y[2])-l_1*u[1])*std::cos(y[2])+(2*A*std::pow(y[1],2)+C*std::pow(y[1],2)*std::cos(y[2])+2*A*y[1]*y[3])*std::sin(y[2])-B*std::sin(y[0]))/(E-C*std::pow(std::cos(y[2]),2));
    yy[2]=y[3];
    yy[3]=(R*u[1]-F*std::sin(y[2])-G*std::sin(y[0])-H*std::cos(y[2])*std::sin(y[2])-Mm*std::cos(y[2])*std::sin(y[0])-Mm*std::sin(y[0]+y[2])-P*std::cos(y[2])-Q*std::sin(y[0]+y[2])*std::cos(y[2])-std::pow(l_2,2)*m_2*u[0])/(S-S*std::pow(std::cos(y[2]),2)+T);

    yy[4]=L[0][1]*y[5];
    yy[5]=L[1][0]*y[4]+L[1][1]*y[5]+L[1][2]*y[6]+L[1][3]*y[7]+w[1];
    yy[6]=L[2][3]*y[7];
    yy[7]=L[3][0]*y[4]+L[3][1]*y[5]+L[3][2]*y[6]+L[3][3]*y[7]+w[3];
  };
  //while(ignore==false){
    scots::runge_kutta_fixed4(rhs,y,u,dis, w3_lb,w3_ub, 2*state_dim,tau,10);
  //  ignore = dis.get_out_of_domain();
};

auto rs_repost = [&dis,&ge,w3_lb,w3_ub](ds_type &y, input_type &u, bool &neigbour) -> void {
  dis.set_intersection_check();
  //dis.set_out_of_domain();
  auto rhs =[&dis,&ge,w3_lb,w3_ub](ds_type &yy, const ds_type &y, input_type &u) -> void {
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
    //disturbance_type w = {0.05, 0.05, 0.05};
    double L[4][4];
   L[0][1]=1;
    L[2][3]=1;
    L[1][0]=6.9296;
    L[1][1]=0.2210;
    L[1][2]=19.7559+2*(u[0]+2*u[1])-1.0005e-07*(u[0]-3*u[1])-2.3020e-07*(u[0]-2*u[1]);
    L[1][3]=0.0628;
    L[3][0]=31.8054;
    L[3][1]=0.4638;
    L[3][2]=9.8780 + 0.1113*(20*u[0]-20*u[1])+0.0369*20*u[1];
    L[3][3]=0.1013;
    double F=(2*l_1*std::pow(l_2,3)*std::pow(m_2,2) + std::pow(l_1,3) *l_2*std::pow(m_2,2) + std::pow(l_1,3) *l_2*m_2*m_1)* std::pow(y[1],2)+2*l_1*std::pow(l_2,3)*std::pow(m_2,2)*y[1]*y[3];
    double H=std::pow(l_1,2)*std::pow(l_2,2)*std::pow(m_2,2)*(3*std::pow(y[2],2)+2*y[1]*y[3]);
    double P=l_1*l_2*m_2*u[0]-2*l_1*l_2*m_2*u[1];
  
    yy[0]=y[1];
    yy[1]=(l_2*(u[0]-u[1])+(D*std::sin(y[0]+y[2])-l_1*u[1])*std::cos(y[2])+(2*A*std::pow(y[1],2)+C*std::pow(y[1],2)*std::cos(y[2])+2*A*y[1]*y[3])*std::sin(y[2])-B*std::sin(y[0]))/(E-C*std::pow(std::cos(y[2]),2));
    yy[2]=y[3];
    yy[3]=(R*u[1]-F*std::sin(y[2])-G*std::sin(y[0])-H*std::cos(y[2])*std::sin(y[2])-Mm*std::cos(y[2])*std::sin(y[0])-Mm*std::sin(y[0]+y[2])-P*std::cos(y[2])-Q*std::sin(y[0]+y[2])*std::cos(y[2])-std::pow(l_2,2)*m_2*u[0])/(S-S*std::pow(std::cos(y[2]),2)+T);

    yy[4]=L[0][1]*y[5];
    yy[5]=L[1][0]*y[4]+L[1][1]*y[5]+L[1][2]*y[6]+L[1][3]*y[7]+w[1];
    yy[6]=L[2][3]*y[7];
    yy[7]=L[3][0]*y[4]+L[3][1]*y[5]+L[3][2]*y[6]+L[3][3]*y[7]+w[3];
   
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
  scots::TransitionFunction tf_o1d;
  scots::TransitionFunction tf_new;
  scots::TransitionFunction tf_standard;
  scots::Abstraction<state_dim,input_dim> abs(ss,is);
  std::queue<abs_type> online_queue;
  //dis.update_disturbance(w_2, w2_lb, w2_ub);
  tt.tic();
  abs.compute_gb(tf_o1d,rs_post);
  //abs.compute_gb(tf,vehicle_post, radius_post);
  tt.toc();

  //if(!getrusage(RUSAGE_SELF, &usage))
 //   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_o1d.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_o1d.get_no_transitions() << std::endl;
    dis.disturbance_readlocal();
   // dis.update_disturbance(w_3, w3_lb, w3_ub);
  std::cout << "\nComputing the new transition function locally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.recompute_gb(tf_new, online_queue, tf_o1d, w3_lb, w3_ub, rs_repost);
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of new transitions: " << tf_new.get_no_transitions() << std::endl;
  double t1=tt.toc();

   std::cout << "\nComputing the stardard transition function globally (after distrubance changes): " << std::endl;
  tt.tic();
  abs.compute_gb(tf_standard,rs_post);
  
 if(!getrusage(RUSAGE_SELF, &usage))
   std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_new.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf_standard.get_no_transitions() << std::endl;
  double t2= tt.toc();
  
  

  /* define target set */
  auto target = [&ss,&s_eta](const abs_type& idx) {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (M_PI/2 <= (x[0]-s_eta[0]/2.0) && (x[0]+s_eta[0]/2.0) <= M_PI/2+0.2 && 
        0 <= (x[2]-s_eta[2]/2.0) && (x[2]+s_eta[2]/2.0) <= 0.2)
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


    // std::ofstream write;
    // std::ifstream read;
    // write.open("result.txt", std::ios::app);          
    // write << "p:"<<persent<<" t1:"<<t1<<" t2:"<<t2<<std::endl;
    // write.close();
    // read.close();

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
 
  for (int k = 24; k < 25; ++k)
  {
    
      main_parameters(k);
  }
  return 1;
}
