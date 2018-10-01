

#ifndef GbEstimation_HH_
#define GbEstimation_HH_

#include <iostream>
#include <cstring>
#include <vector>
#include <queue>

#include "UniformGrid.hh"


namespace scots{

template<class disturbance_type>
class GbEstimation
{
private:
	disturbance_type w_max;
	disturbance_type w_integrals;
	std::vector<disturbance_type> l_exp;
	int ss_dim;
	abs_type M;
public:
	GbEstimation(const UniformGrid& is, const UniformGrid& ss,const disturbance_type w1, const disturbance_type w2){
		ss_dim=ss.get_dim();
		M = is.size();
		for (int i = 0; i < ss_dim; ++i)
		{
			w_max[i]=std::max(w1[i],w2[i]);
			w_integrals[i]=0;
		}
		
		l_exp.resize(M);
	}
	~GbEstimation();

	template<class F>
	void exp_interals(F& l_matrix, const double t){
	
		auto f=[M,t,l_matrix,ss_dim,w_max](double x){
		disturbance_type le = scots::matrixExp(M,l_matrix,t-x);
		disturbance_type f;
		for (int i = 0; i < ss_dim; ++i)
			f[i]=0;
			for(int j=0; j<ss_dim; ++j)
				f[i]+=le[i][j]*w_max[j];
		
		return f;
		};

		w_integrals=scots::simpson3(f,0,t, 6, ss_dim);
		l_exp=scots::matrixExp(M,l_matrix,t);
	}

	template<class F1, class F2>
	disturbance_type gb_estimate(F1& r, F2& u){
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
