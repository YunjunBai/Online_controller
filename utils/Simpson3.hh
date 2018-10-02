

#ifndef SIMPSON3_HH_
#define SIMPSON3_HH_

#include "MatrixExp.hh"
namespace scots{

template<class F1, class disturbance_type, class F3>
void simpson3(F1& f, disturbance_type& integral, F3 input_id, const double a, const double b, const int n, const int dim) noexcept {
    
    double h=(b-a)/n; 
    disturbance_type sum;

    
      double x[n+1];
      disturbance_type y[n+1];                     
      for (int i=0;i<n+1;i++)
      {                        //loop to evaluate x0,...xn and y0,...yn
          x[i]=a+i*h;                //and store them in arrays
          y[i]=f(x[i],input_id);
      }
      for (int i=1;i<n;i+=2)
      {
        for(int j=0;j<dim;j++)
          sum[j]=sum[j]+4.0*y[i][j];
                          //loop to evaluate 4*(y1+y3+y5+...+yn-1)
      }
      for (int i=2;i<n-1;i+=2)
      {
        for (int j = 0; j < dim; j++)
          sum[j]=sum[j]+2.0*y[i][j];                /*loop to evaluate 4*(y1+y3+y5+...+yn-1)+
                                             2*(y2+y4+y6+...+yn-2)*/ 
      } 
      for (int j = 0; j < dim; j++)
        integral[j]=h/3.0*(y[0][j]+y[n][j]+sum[j]);    //h/3*[y0+yn+4*(y1+y3+y5+...+yn-1)+2*(y2+y4+y6+...+yn-2)]
    

}

}
#endif