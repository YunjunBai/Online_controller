

#ifndef SIMPSON3_HH_
#define SIMPSON3_HH_

namespace scots{

template<class F, class ds_type>
disturbance_type simpson3(F& f,double& a, double& b, int& n, int& dim) noexcept {
   
    ds_type integral;
    
    double h=(b-a)/n;  
    for (int j = 0; j < dim; ++j){
      double x[n+1],y[n+1];                     /get the width of the subintervals
      for (i=0;i<n+1;i++)
      {                        //loop to evaluate x0,...xn and y0,...yn
          x[i]=a+i*h;                //and store them in arrays
          y[i]=f(x[i]);
      }
      for (i=1;i<n;i+=2)
      {
          sum=sum+4.0*y[i];                //loop to evaluate 4*(y1+y3+y5+...+yn-1)
      }
      for (i=2;i<n-1;i+=2)
      {
          sum=sum+2.0*y[i];                /*loop to evaluate 4*(y1+y3+y5+...+yn-1)+
                                             2*(y2+y4+y6+...+yn-2)*/ 
      } 
      integral[j]=h/3.0*(y[0]+y[n]+sum);    //h/3*[y0+yn+4*(y1+y3+y5+...+yn-1)+2*(y2+y4+y6+...+yn-2)]
    }
    return integral;
}

}
#endif