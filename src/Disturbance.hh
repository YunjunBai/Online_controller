/*disturbance class*/

#ifndef DISTURBANCE_HH_
#define DISTURBANCE_HH_

#include <iostream>
#include <cstring>
#include <vector>
#include <queue>

#include "UniformGrid.hh"


namespace scots{

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

public:
	/*initial constructor, w is a subset of W, initially for all states in X , the disturbance is w*/
	Disturbance( const disturbance_type w,
				 const UniformGrid& ss):disturbance_(w), states_alphabet(ss){
		abs_type N = states_alphabet.size();
		x_disturbance.reset(new disturbance_type[N]);
		for (abs_type i=0; i< N; ++i)
			x_disturbance[i] = disturbance_;	
		init_disturbance=disturbance_;
	}

	~Disturbance()=default;

	/*BDD subset_X is the region where the disturbance need to change to new_w*/
	template<class F1, class F2>
	void update_disturbance(F1& new_disturbance, 
							F2& d_lb, 
							F2& d_ub){
		/*get all indices in that region, put ids into a queue*/
		std::queue<abs_type> id_queue;
		int ss_dim = states_alphabet.get_dim();
		std::vector<abs_type> NN=states_alphabet.get_nn();
		std::vector<double>  ss_ll=states_alphabet.get_lower_left();
		std::vector<double> ss_ur=states_alphabet.get_upper_right();
		std::vector<double> eta = states_alphabet.get_eta();
		std::vector<abs_type> lb(ss_dim);  /* lower-left corner */
    	std::vector<abs_type> ub(ss_dim);  /* upper-right corner */
    	std::vector<abs_type> no(ss_dim);  /* number of cells per dim */
    	std::vector<abs_type> cc(ss_dim);  /* coordinate of current cell in the post */

		abs_type nRegion=1;
		for (int i = 0; i < ss_dim; ++i)
		{
          /* check for out of bounds */
          double left = d_lb[i];
          double right = d_ub[i];
          if(left <= ss_ll[i]-eta[i]/2.0  || right >= ss_ur[i]+eta[i]/2.0)  {
            //out_of_domain[i*M+j]=true;
            break;
          } 

          /* integer coordinate of lower left corner of post */
          lb[i] = static_cast<abs_type>((left-ss_ll[i]+eta[i]/2.0)/eta[i]);
          /* integer coordinate of upper right corner of post */
          ub[i] = static_cast<abs_type>((right-ss_ll[i]+eta[i]/2.0)/eta[i]);
          /* number of grid points in the post in each dimension */
          no[i]=(ub[i]-lb[i]+1);
          /* total number of post */
          nRegion*=no[i];
          cc[i]=0;
        }

        /* compute indices */
        for(abs_type k=0; k<nRegion; k++) {
          abs_type q=0;
          for(int l=0; l<ss_dim; l++) 
            q+=(lb[l]+cc[l])*NN[l];
          cc[0]++;
          for(int l=0; l<ss_dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
          id_queue.push(q);
		}
		/*trasfer all ids of subss to change the values in x_disturbance vector */
		while(!id_queue.empty()){
			abs_type q=id_queue.front();
			id_queue.pop();
			x_disturbance[q]=new_disturbance;
		}
	}

	/*given a x, return the related w*/
	template<class F1, class F2>
    disturbance_type get_disturbance(F1& x, F2& r){

        int ss_dim = states_alphabet.get_dim();
        std::vector<abs_type> m_NN = states_alphabet.get_nn();
 
        disturbance_type max_w=init_disturbance;
        std::vector<abs_type> lb(ss_dim);  /* lower-left corner */
    	std::vector<abs_type> ub(ss_dim);  /* upper-right corner */
    	std::vector<abs_type> no(ss_dim);  /* number of cells per dim */
    	std::vector<abs_type> cc(ss_dim);  /* coordinate of current cell in the (x,r) */
        state_type lower_left;
    	state_type upper_right;
    	state_type m_eta;
    	state_type m_z;
        for(int i=0; i<ss_dim; i++) {
      		m_eta[i]=states_alphabet.get_eta()[i];
      		m_z[i]=states_alphabet.get_eta()[i]/1e10;
      		lower_left[i]=states_alphabet.get_lower_left()[i];
      		upper_right[i]=states_alphabet.get_upper_right()[i];
      	}
        abs_type nintersection=1;
        for(int k=0; k<ss_dim; k++) {
          /* check for out of bounds */
          double left = x[k]-r[k]-m_z[k];
          double right = x[k]+r[k]+m_z[k];
          if(left <= lower_left[k]-m_eta[k]/2.0  || right >= upper_right[k]+m_eta[k]/2.0)  {
            return max_w;
          } 

          /* integer coordinate of lower left corner of (x,r) */
          lb[k] = static_cast<abs_type>((left-lower_left[k]+m_eta[k]/2.0)/m_eta[k]);
          /* integer coordinate of upper right corner of (x,r) */
          ub[k] = static_cast<abs_type>((right-lower_left[k]+m_eta[k]/2.0)/m_eta[k]);
          /* number of grid points in the post in each dimension */
          no[k]=(ub[k]-lb[k]+1);
          /* total number of (x,r) */
          nintersection*=no[k];
          cc[k]=0;
        }
        /* compute indices */
        for(abs_type k=0; k<nintersection; k++) {
          abs_type q=0;
          for(int l=0; l<ss_dim; l++){
          	q+=(lb[l]+cc[l])*m_NN[l];
          	max_w[l] = (x_disturbance[q][l] > max_w[l] ? x_disturbance[q][l] : max_w[l]);
          } 
            
          cc[0]++;
          for(int l=0; l<ss_dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
        }
    return max_w;
    		
	}
	
};//class closed
}//namespace closed

#endif