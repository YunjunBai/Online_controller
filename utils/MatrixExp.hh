

#ifndef MATRIXEXP_HH_
#define MATRIXEXP_HH_

namespace scots{

template<class F>
disturbance_type matrixExp(abs_type M, F& l_matrix){


	for (int i = 0; i < M; ++i)
	{
		double l[ss_dim][ss_dim]=l_matrix(i);

	}


	return ;
}

}
#endif