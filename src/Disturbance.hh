/*disturbance class*/

#ifndef DISTURBANCE_HH_
#define DISTURBANCE_HH_

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <queue>

#include "UniformGrid.hh"

namespace scots{
namespace params {
  auto avoid_dis = [](const abs_type&) noexcept {return false;};
}

template<class disturbance_type, class state_type>
class Disturbance
{
private:
  /*real values*/
  disturbance_type disturbance_;
  UniformGrid states_alphabet;
  disturbance_type init_disturbance;
  /*vector containing the disturbance for each grid state i*/
  std::unique_ptr<disturbance_type[]>  x_disturbance;
  std::vector<abs_type> pre_CornerIDs;
  disturbance_type pre_max_w;
  /*mark if disturbance already updates or not*/
  bool disturbance_marker;
  int ss_dim;
  std::vector<abs_type> no;  /* number of cells per dim */
  std::vector<abs_type> cc;  /* coordinate of current cell in the (x,r) */
  std::vector<abs_type> idx; 
  std::vector<abs_type> base;
  std::vector<abs_type> m_NN; 
  std::vector<double> lower_left; 
  std::vector<double> upper_right; 
  std::vector<double> m_eta;

  bool intersection_check;
  bool out_of_domain;

public:
  /*initial constructor, w is a subset of W, initially for all states in X , the disturbance is w*/
  Disturbance( const disturbance_type w,
         const UniformGrid& ss):disturbance_(w), states_alphabet(ss){
    abs_type N = states_alphabet.size();
    x_disturbance.reset(new disturbance_type[N]);
    for (abs_type i=0; i< N; ++i)
      x_disturbance[i] = disturbance_;  
    init_disturbance=disturbance_;
    pre_CornerIDs={2*N,N};
    pre_max_w=init_disturbance;
    disturbance_marker=false;
    ss_dim = states_alphabet.get_dim();
    no.resize(ss_dim);  /* number of cells per dim */
    cc.resize(ss_dim);  /* coordinate of current cell in the (x,r) */
    idx.resize(ss_dim); 
    base.resize(ss_dim);
    m_NN = states_alphabet.get_nn();
    lower_left = states_alphabet.get_lower_left();
    upper_right = states_alphabet.get_upper_right();
    m_eta = states_alphabet.get_eta();
    intersection_check=false;
    out_of_domain=false;
  }

  ~Disturbance()=default;

  void disturbance_read(){
    std::string line;
    std::ifstream myfile("DistMatGlobal_1.txt");
    double value[ss_dim*3];
    disturbance_type d;
    disturbance_type lb;
    disturbance_type ub;
  for (std::string line; std::getline(myfile,line);  ){
    std::string::size_type sz;    
    for (int i = 0; i < 12; ++i)
    {
      double tmp=std::stod(line,&sz);
      value[i]=tmp;
    }
   for (int j = 0; j < ss_dim; ++j)
      {
        d[j]=value[j];
        lb[j]=value[j+ss_dim];
        ub[j]=value[j+ss_dim*2];
      }
      update_disturbance(d,lb,ub);
      //value.clear();
  }
  }

  void disturbance_readlocal(){
    std::string line;
    std::ifstream myfile("DistMatLocal_1.txt");
    double value[ss_dim*3];
    disturbance_type d;
    disturbance_type lb;
    disturbance_type ub;
  for (std::string line; std::getline(myfile,line);  ){
    std::string::size_type sz;    
    for (int i = 0; i < 12; ++i)
    {
      double tmp=std::stod(line,&sz);
      value[i]=tmp;
    }
   for (int j = 0; j < ss_dim; ++j)
      {
        d[j]=value[j];
        lb[j]=value[j+ss_dim];
        ub[j]=value[j+ss_dim*2];
      }
      manipulator_local(d,lb,ub);
      //value.clear();
  }
  }
  /*BDD subset_X is the region where the disturbance need to change to new_w*/
  template<class F1, class F2, class F4=decltype(params::avoid_dis)>
  void update_disturbance(F1& new_disturbance, 
              F2& d_lb, 
              F2& d_ub,
              F4& avoid=params::avoid_dis){

    disturbance_marker=true;
    /*get all indices in that region, put ids into a queue*/
    //std::queue<abs_type> id_queue;
    std::vector<abs_type> lb(ss_dim,0);  /* lower-left corner */
    std::vector<abs_type> ub(ss_dim,0);  /* upper-right corner */
    std::vector<abs_type> no(ss_dim,0);  /* number of cells per dim */
    std::vector<abs_type> cc(ss_dim,0);  /* coordinate of current cell in the post */
     //std::unique_ptr<bool[]> recomputed_mark(new bool[N*M]());
    abs_type nRegion=1;
    for (int i = 0; i < ss_dim; ++i)
    {
          /* check for out of bounds */
          double left = d_lb[i]-m_eta[i]/1e10;
          double right = d_ub[i]+m_eta[i]/1e10;
          if(left <= lower_left[i]-m_eta[i]/2.0)
            left=lower_left[i];
          if (right >= upper_right[i]+m_eta[i]/2.0)
            right=upper_right[i];

          /* integer coordinate of lower left corner of post */
          lb[i] = static_cast<abs_type>((left-lower_left[i]+m_eta[i]/2.0)/m_eta[i]);
          /* integer coordinate of upper right corner of post */
          ub[i] = static_cast<abs_type>((right-lower_left[i]+m_eta[i]/2.0)/m_eta[i]);
          /* number of grid points in the post in each dimension */
          no[i]=(ub[i]-lb[i]+1);
          /* total number of post */
          nRegion*=no[i];
          cc[i]=0;
        }
       // std::cout<<"The number of states in the new disturbance region: "<<nRegion<<std::endl;
        /* compute indices */
        for(abs_type k=0; k<nRegion; k++) {
          abs_type q=0;
          for(int l=0; l<ss_dim; l++) 
            q+=(lb[l]+cc[l])*m_NN[l];
          cc[0]++;
          for(int l=0; l<ss_dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
          if(avoid(q))
            continue;
         x_disturbance[q]=new_disturbance;
    }
  
  }
  disturbance_type get_initDisturbance(){
    return init_disturbance;
  }

  /*given a x, return the related w*/
  template<class F1, class F2, class F4=decltype(params::avoid_dis)>
    disturbance_type get_disturbance(F1& x, F2& r, F4& avoid=params::avoid_dis){
        abs_type no[ss_dim];
        abs_type cc[ss_dim];
        abs_type idx[ss_dim];
        abs_type base[ss_dim];
        /*first check if disturbance update or not, if no, return initial disturbance*/
        if(!disturbance_marker){
          return init_disturbance;
        }

        std::vector<abs_type> corner_IDs(2,0);
        // for(int k=0; k<ss_dim; k++){
        //   no[k]=0;
        //   base[k]=0;
        // }
        for(int k=0; k<ss_dim; k++) {
            /* check for out of bounds */
            double left = x[k]-r[k]-m_eta[k]/1e10;
            double right = x[k]+r[k]+m_eta[k]/1e10;
            if (left <= lower_left[k]-m_eta[k]/2.0 || left>=upper_right[k]+m_eta[k]/2.0 
                || right >= upper_right[k]+m_eta[k]/2.0 || right <= lower_left[k]-m_eta[k]/2.0){
              out_of_domain = true;
              return init_disturbance;
              break;
            }

            /* integer coordinate of lower left corner of (x,r) */
            abs_type lb = static_cast<abs_type>((left-lower_left[k]+m_eta[k]/2.0)/m_eta[k]);
            /* integer coordinate of upper right corner of (x,r) */
            abs_type ub = static_cast<abs_type>((right-lower_left[k]+m_eta[k]/2.0)/m_eta[k]);

           /* number of grid points in the post in each dimension */
            no[k] = (ub-lb+1);
            
            base[k] = lb * m_NN[k];
            corner_IDs[0] += lb * m_NN[k]; //low-left id
            corner_IDs[1] += ub * m_NN[k]; //up-rigth id
            cc[k] = 0;
            idx[k] = 0;
        }

        if(pre_CornerIDs[0]==corner_IDs[0] && pre_CornerIDs[1]==corner_IDs[1])
          return pre_max_w; 
        else{
          pre_CornerIDs[0]=corner_IDs[0];
          pre_CornerIDs[1]=corner_IDs[1];
        
          /* compute indices */
          int i = ss_dim - 1;
          cc[i] = base[i];
          /*initialize max_w with the disturbance of first grid piont in this region*/
          disturbance_type max_w;
          max_w=x_disturbance[corner_IDs[0]];
          
          while(i < ss_dim) {
              if (idx[i] < no[i]) {
                  if (i > 0) { //not yet computed the index for all the dimensions
                      cc[i-1] = cc[i] + base[i-1];
                      i -= 1;
                  } else { //we have all the dimensions except 0, we can look the disturbance up
                      abs_type q = cc[0];
                      for (abs_type j = 0; j < no[0]; j++) {
                          if(avoid(q))
                            continue;
                          for(int l=0; l<ss_dim; l++){
                              max_w[l] = std::max(x_disturbance[q][l], max_w[l]);
                          }
                          q += m_NN[0];
                      }
                      idx[0] = no[0]; // make it backtrack to
                  }
              } else {
                  //we reached the bound, going back to the previous dimension
                  idx[i] = 0;
                  i += 1;
                  if (i < ss_dim) {
                      idx[i] += 1;
                      cc[i] += m_NN[i];
                  }
              }
          }
          pre_max_w = max_w;
          
          return max_w;
        }
  }
  
template<class ds_type>
void intersection(ds_type y, state_type d_lb,state_type d_ub){
  bool tmp=true;
  state_type x;
  state_type r;
  for (int k = 0; k < ss_dim; k++)
  {
    x[k]=y[k];
    r[k]=y[k+ss_dim];
    /* integer coordinate of lower left corner of (x,r) */
  abs_type lb = static_cast<abs_type>((d_lb[k]-lower_left[k]+m_eta[k]/2.0)/m_eta[k]);
  /* integer coordinate of upper right corner of (x,r) */
  abs_type ub = static_cast<abs_type>((d_ub[k]-lower_left[k]+m_eta[k]/2.0)/m_eta[k]);
  d_lb[k]=m_eta[k] * lb + lower_left[k];
  d_ub[k]=m_eta[k] * (ub +1) + lower_left[k];
  }
  
  for(int i=0; i<ss_dim; i++){
    if ((x[i]-r[i])>d_ub[i]
      || d_lb[i] >(x[i]+r[i]))
    {
      tmp=false;
      break;
    }
  }
  if(tmp)
    intersection_check = true;  
}

state_type get_maxdistance(state_type max_dynamic, const double tau){
  state_type max_distance;
  for (int i = 0; i < ss_dim; ++i)
  {
    max_distance[i]=max_dynamic[i]*tau+m_eta[i];
  }
  return max_distance;
}

bool get_intersection_check(){
  return intersection_check;
}
 void set_intersection_check(){
  intersection_check=false;
 }
bool get_out_of_domain(){
  return out_of_domain;
}


void set_out_of_domain(){
  out_of_domain=false;
}

template<class F1, class F2, class F4=decltype(params::avoid_dis)>
  void manipulator_local(F1& new_disturbance, 
              F2& d_lb, 
              F2& d_ub,
              F4& avoid=params::avoid_dis){
     abs_type N = states_alphabet.size();
    disturbance_marker=true;
    /*get all indices in that region, put ids into a queue*/
    //std::queue<abs_type> id_queue;
    std::vector<abs_type> lb(ss_dim,0);  /* lower-left corner */
    std::vector<abs_type> ub(ss_dim,0);  /* upper-right corner */
    std::vector<abs_type> no(ss_dim,0);  /* number of cells per dim */
    std::vector<abs_type> cc(ss_dim,0);  /* coordinate of current cell in the post */
    std::unique_ptr<bool[]> update_marker(new bool[N]());
    abs_type nRegion=1;
    for (int i = 0; i < ss_dim; ++i)
    {
          /* check for out of bounds */
          double left = d_lb[i]-m_eta[i]/1e10;
          double right = d_ub[i]+m_eta[i]/1e10;
          if(left <= lower_left[i]-m_eta[i]/2.0)
            left=lower_left[i];
          if (right >= upper_right[i]+m_eta[i]/2.0)
            right=upper_right[i];

          /* integer coordinate of lower left corner of post */
          lb[i] = static_cast<abs_type>((left-lower_left[i]+m_eta[i]/2.0)/m_eta[i]);
          /* integer coordinate of upper right corner of post */
          ub[i] = static_cast<abs_type>((right-lower_left[i]+m_eta[i]/2.0)/m_eta[i]);
          /* number of grid points in the post in each dimension */
          no[i]=(ub[i]-lb[i]+1);
          /* total number of post */
          nRegion*=no[i];
          cc[i]=0;
        }
       // std::cout<<"The number of states in the new disturbance region: "<<nRegion<<std::endl;
        /* compute indices */
        for(abs_type k=0; k<nRegion; k++) {
          abs_type q=0;
          for(int l=0; l<ss_dim; l++) 
            q+=(lb[l]+cc[l])*m_NN[l];
          cc[0]++;
          for(int l=0; l<ss_dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
          if(avoid(q))
            continue;
          if (!update_marker[q])
          {
            x_disturbance[q]=new_disturbance;
          }
          else{
            x_disturbance[q]=std::max(x_disturbance[q],new_disturbance);
          }
         
    }
  }
	
};//class closed
}//namespace closed

#endif
