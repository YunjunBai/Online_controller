

#ifndef MATRIXEXP_HH_
#define MATRIXEXP_HH_
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;
namespace scots{

template<class F1,  class matrix_type>
void matrixExp(F1& l_matrix, matrix_type& l_exp, const abs_type M, abs_type input_id, const double t, const int ss_dim){

		MatrixXd l;
		l=l_matrix(input_id);
		for (int j = 0; j < ss_dim; ++j)
			for (int k = 0; k < ss_dim; ++k)
				l(j,k)*=t;

		l_exp=l.exp();
	
}

}
#endif