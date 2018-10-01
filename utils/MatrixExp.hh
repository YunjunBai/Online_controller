

#ifndef MATRIXEXP_HH_
#define MATRIXEXP_HH_
#include <unsupported/Eigen/MatrixFunctions>
namespace scots{

template<class F>
disturbance_type matrixExp(abs_type M, F& l_matrix, double& t){

	std::vector<disturbance_type> l_exp;
	for (int i = 0; i < M; ++i)
	{
		double l[ss_dim][ss_dim]=l_matrix(i);
		for (int j = 0; j < ss_dim; ++j)
			for (int k = 0; k < ss_dim; ++k)
				l[j][k]*=t;

		l_exp[i]=l.exp();
	}

	return l_exp;
}

}
#endif