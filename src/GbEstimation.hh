

#ifndef GbEstimation_HH_
#define GbEstimation_HH_

#include <iostream>
#include <cstring>
#include <vector>
#include <queue>
#include "Simpson3.hh"
#include "MatrixExp.hh"
#include "UniformGrid.hh"


namespace scots{

template<class disturbance_type, class matrix_type>
class GbEstimation
{
private:
	disturbance_type w_max;
	disturbance_type w_integrals;
	std::vector<matrix_type> l_exp;
	int ss_dim;
	abs_type M;
	UniformGrid input_space;
public:
	GbEstimation(const UniformGrid& is, const UniformGrid& ss,const disturbance_type w1, const disturbance_type w2){
		ss_dim=ss.get_dim();
		M = is.size();
		input_space=is;
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
		const abs_type M_=M;
		int ss_dim_=ss_dim;
		disturbance_type w_max_=w_max;

		auto f=[M_,ss_dim_,w_max_,t,l_matrix](const double x,  abs_type input_id){
			matrix_type le;
			scots::matrixExp(l_matrix, le, M_, input_id, t-x, ss_dim_);
			disturbance_type f_value;
			for (int i = 0; i < ss_dim_; ++i){
				f_value[i]=0;
				for(int j=0; j<ss_dim_; ++j)
					f_value[i]+=le[i][j] * w_max_[j];
			}
			return f_value;
		};
		for (abs_type i = 0; i < M; ++i)
		{
			scots::simpson3(f,w_integrals,i,0,t, 6, ss_dim);
			scots::matrixExp(l_matrix, l_exp[i], M, i, t, ss_dim);
		}
		
	}

	template<class F1, class F2>
	disturbance_type gb_estimate(F1& r, F2& u){
		disturbance_type r_es;
		abs_type u_id=input_space.xtoi(u);
		for (int i = 0; i < ss_dim; i++){
			r_es[i]=0;
			for (int j = 0; j < ss_dim; j++)
				r_es[i]=l_exp[u_id][i][j]*r[i]+w_integrals[i];
			
		}
			
		
		return r_es;
	}
};

}//namespace closed

#endif
