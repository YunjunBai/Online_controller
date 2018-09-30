

#ifndef GbEstimation_HH_
#define GbEstimation_HH_

#include <iostream>
#include <cstring>
#include <vector>
#include <queue>

#include "UniformGrid.hh"


namespace scots{

template<class disturbance_type, class input_type>
class GbEstimation
{
private:
	disturbance_type w_max;
	disturbance_type w_integrals;
	std::vector<disturbance_type> l_exp;
public:
	GbEstimation(const UniformGrid& is,const UniformGrid& ss, const disturbance_type w1, const disturbance_type w2){
		int ss_dim=ss.get_dim();
		abs_type M = is.size();
		for (int i = 0; i < ss_dim; ++i)
		{
			w_max[i]=std::max(w1[i],w2[i]);
		}

		w_integrals=scots::simpson3(,w_max);
		l_exp=scots::matrixExp(M,l_matrix);
	}

	~GbEstimation();

	disturbance_type gb_estimate(state_type r, input_type u){
		state_type r_es;
		abs_type u_id=is.xtoi(u);
		for (int i = 0; i < ss_dim; ++i)
		{
			r_es[i]=a[u_id][i]*r[i]+b[i];
		}
		return r_es;
	}

};

}//namespace closed

#endif
