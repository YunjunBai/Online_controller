

#ifndef MATRIXEXP_HH_
#define MATRIXEXP_HH_
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;
namespace scots{

template<class F1,  class matrix_type>
void matrixExp(F1& l_matrix, matrix_type& l_exp, abs_type input_id, const double t, const int ss_dim){

		matrix_type l;
		MatrixXd l_(ss_dim,ss_dim);
		MatrixXd l_exp_;
		l=l_matrix(input_id);
		for (int j = 0; j < ss_dim; ++j)
			for (int k = 0; k < ss_dim; ++k){
				l[j][k]*=t;
				l_(j,k)=l[j][k];
			}

		l_exp_=l_.exp();
		for (int j = 0; j < ss_dim; ++j)
			for (int k = 0; k < ss_dim; ++k)
				l_exp[j][k]=l_exp_(j,k);
	
}

}
#endif