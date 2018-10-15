/*
 * Abstraction.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/** @file **/
#ifndef ABSTRACTION_HH_
#define ABSTRACTION_HH_

#include <iostream>
#include <cstring>
#include <memory>
#include <vector>
#include "UniformGrid.hh"
#include "TransitionFunction.hh"
#include "Disturbance.hh"
#include "InputOutput.hh"

/** @namespace scots **/ 
namespace scots {

/** @cond **/
/* default parameter for the third parameter of Abstraction::compute_gb  */
namespace params {
  auto avoid_abs = [](const abs_type&) noexcept {return false;};
}
/** @endcond **/


/**
 * @class Abstraction
 * 
 * @brief Computation of the transition function of symbolic model/abstraction
 *
 * It implements in compute_gb the computation of the transition function based
 * on a growth bound. Additionally, it provides the functions print_post and
 * print_post_gb to print the center of the cells that cover the attainable set
 * which is computed from the growth bound.
 *
 * See 
 * - the manual in <a href="./../../manual/manual.pdf">manual.pdf</a>
 * - http://arxiv.org/abs/1503.03715 for theoretical background 
 *
 **/
template<class state_type, class input_type, class ds_type>
class Abstraction:public TransitionFunction {
private:
  /* grid information of state alphabet */
  const UniformGrid m_state_alphabet;
  /* grid information of input alphabet */
  const UniformGrid m_input_alphabet;
  /* measurement error bound */
  std::unique_ptr<double[]> m_z;
  /* print progress to the console (default m_verbose=true) */
  bool m_verbose=true;
  /* to display the progress of the computation of the abstraction */
  void progress(const abs_type& i, const abs_type& N, abs_type& counter) {
    if(!m_verbose)
      return;
    if(((double)i/(double)N*100)>=counter){
      if((counter%10)==0)
        std::cout << counter;
      else if((counter%2)==0) 
        std::cout << ".";
      counter++;
    }
    std::flush(std::cout); 
    if(i==(N-1))
      std::cout <<"100\n";
    
  }
  void progress_re(const abs_type& N,const abs_type& M, abs_type& counter, abs_type& counter_states) {
    if(!m_verbose)
      return;
    if(((double)counter_states/(double)N*M*100)>=counter){
      if((counter%10)==0)
        std::cout << counter;
      else if((counter%2)==0) 
        std::cout << ".";
      counter++;
    }
    counter_states++;
    std::flush(std::cout);    
  }

public:
  /* @cond  EXCLUDE from doxygen*/
  /* deactivate standard constructor */
  Abstraction() = delete;
  /* cannot be copied or moved */
  Abstraction(Abstraction&&) = delete;
  Abstraction(const Abstraction&) = delete;
  Abstraction& operator=(Abstraction&&)=delete;
  Abstraction& operator=(const Abstraction&)=delete;
  /* @endcond */

  /** @brief constructor with the abstract state alphabet and abstract input
   * alphabet (measurement error bound is set to zero)
   *
   *  @param state_alphabet - UniformGrid representing the abstract state alphabet
   *  @param input_alphabet - UniformGrid representing the abstract input alphabet
   **/
  Abstraction(const UniformGrid& state_alphabet,
              const UniformGrid& input_alphabet) :
              m_state_alphabet(state_alphabet),
              m_input_alphabet(input_alphabet),
              m_z(new double[state_alphabet.get_dim()]()) {
  /* default value of the measurement error 
     * (heurisitc to prevent rounding errors)*/
    for(int i=0; i<m_state_alphabet.get_dim(); i++)
      m_z[i]=m_state_alphabet.get_eta()[i]/1e10;
  }

  /** 
   * @brief computes the transition function
   *
   * @param[out] transition_function - the result of the computation
   *
   * @param[in] system_post - lambda expression of the form
   *                          \verbatim [] (state_type &x, const input_type &u) ->  void  \endverbatim
   *                          system_post(x,u) provides a numerical approximation of ODE 
   *                          solution at time tau with initial state x and input u \n
   *                          the result is stored in x
   *
   * @param[in] radius_post - lambda expression of the form
   *                          \verbatim [] (state_type &r, const state_type& x, const input_type &u) -> void  \endverbatim
   *                          radius_post(x,u) provides a numerical approximation of
   *                          the growth bound for the cell (with center x, radius  r) and input u\n
   *                          the result is stored in r
   *
   * @param[in] avoid  - OPTIONALLY provide lambda expression of the form
   *                     \verbatim [] (const abs_type &i) -> bool \endverbatim
   *                     returns true if the abstract state i is in the avoid
   *                     set; otherwise returns false
   *
   * The computation proceeds in two loops. In the first loop the cornerIDs are
   * comuted, which represent the cell IDs that cover the over-approximation of the attainable set:
   * \verbatim corner_IDs[i*M+j+0] \endverbatim = lower-left cell ID of the
   * integer hyper-interval of cell IDs that cover the over-approximation of the
   * attainable set of cell with ID=i under input ID=j
   * \verbatim corner_IDs[i*M+j+1] \endverbatim = upper-right  of the
   * integer hyper-interval of cell IDs that cover the over-approximation of the
   * attainable set of cell with ID=i under input ID=j
   * 
   * In the second loop the data members of the TransitionFunction are computed.
   * 
   **/
  template<class F1,  class F4=decltype(params::avoid_abs)>
  void compute_gb(TransitionFunction& transition_function, 
                  F1& rs_post, 
                  F4& avoid=params::avoid_abs) {
   
    /* number of cells */
    abs_type N=m_state_alphabet.size(); 
    /* number of inputs */
    abs_type M=m_input_alphabet.size();
    /* number of transitions (to be computed) */
    abs_ptr_type T=0; 
    /* state space dimension */
    int dim=m_state_alphabet.get_dim();
    /* for display purpose */
    abs_type counter=0;
    /* some grid information */
    std::vector<abs_type> NN=m_state_alphabet.get_nn();
    std::queue<abs_type> compute_queue; 
    /* variables for managing the post */
    std::vector<abs_type> lb(dim);  /* lower-left corner */
    std::vector<abs_type> ub(dim);  /* upper-right corner */
    std::vector<abs_type> no(dim);  /* number of cells per dim */
    std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
    /* radius of hyper interval containing the attainable set */
    state_type eta;
    state_type r;
    /* state and input variables */
    state_type x;
    input_type u;
    ds_type y;
    /* for out of bounds check */
    state_type lower_left;
    state_type upper_right;
    /* copy data from m_state_alphabet */
    for(int i=0; i<dim; i++) {
      eta[i]=m_state_alphabet.get_eta()[i];
      lower_left[i]=m_state_alphabet.get_lower_left()[i];
      upper_right[i]=m_state_alphabet.get_upper_right()[i];
    }
    /* init in transition_function the members no_pre, no_post, pre_ptr */ 
    transition_function.init_infrastructure(N,M);
    /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
   
    /*
     * first loop: compute corner_IDs:
     * corner_IDs[i*M+j][0] = lower-left cell index of over-approximation of attainable set 
     * corner_IDs[i*M+j][1] = upper-right cell index of over-approximation of attainable set 
     */
    /* loop over all cells */
    for(abs_type i=0; i<N; i++) {
      /* is i an element of the avoid symbols ? */
      if(avoid(i)) {
        for(abs_type j=0; j<M; j++) {
          transition_function.out_of_domain[i*M+j]=true;
        }
        continue;
      }
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++) {
        transition_function.out_of_domain[i*M+j]=false;
        /* get center x of cell */
        m_state_alphabet.itox(i,x);
        /* cell radius (including measurement errors) */
        for(int k=0; k<dim; k++){
          r[k]=eta[k]/2.0+m_z[k];
          y[k]=x[k];
          y[k+dim]=r[k];
        }
        
        /* current input */
        m_input_alphabet.itox(j,u);
        /* integrate system and radius growth bound */
        /* the result is stored in x and r */
 
      
        rs_post(y,u);
        
        for (int k = 0; k<dim; ++k)
        {
          x[k]=y[k];
          r[k]=y[k+dim];
        }
        
        
        /* determine the cells which intersect with the attainable set: 
         * discrete hyper interval of cell indices 
         * [lb[0]; ub[0]] x .... x [lb[dim-1]; ub[dim-1]]
         * covers attainable set 
         */
        
        abs_type npost=1;
        for(int k=0; k<dim; k++) {
          /* check for out of bounds */
          double left = x[k]-r[k]-m_z[k];
          double right = x[k]+r[k]+m_z[k];
          if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0)  {
            transition_function.out_of_domain[i*M+j]=true;
            break;
          } 

          /* integer coordinate of lower left corner of post */
          lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
          /* integer coordinate of upper right corner of post */
          ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
          /* number of grid points in the post in each dimension */
          no[k]=(ub[k]-lb[k]+1);
          /* total number of post */
          npost*=no[k];
          cc[k]=0;
        }
        //transition_function.corner_IDs[i*(2*M)+2*j]=0;
        //transition_function.corner_IDs[i*(2*M)+2*j+1]=0;
        if(transition_function.out_of_domain[i*M+j]) 
          continue;

        /* compute indices of post */
        for(abs_type k=0; k<npost; k++) {
          abs_type q=0;
          for(int l=0; l<dim; l++) 
            q+=(lb[l]+cc[l])*NN[l];
          cc[0]++;
          for(int l=0; l<dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
          /* (i,j,q) is a transition */    
          /* increment number of pres for (q,j) */ 
          transition_function.m_no_pre[q*M+j]++;
          /* store id's of lower-left and upper-right cell */
          if(k==0)
            transition_function.corner_IDs[i*(2*M)+2*j]=q;
          if(k==npost-1)
            transition_function.corner_IDs[i*(2*M)+2*j+1]=q;
        }
        /* increment number of transitions by number of post */
        T+=npost;
        transition_function.m_no_post[i*M+j]=npost;
      }
      /* print progress */
      if(m_verbose) {
        if(counter==0)
          std::cout << "1st loop: ";
      }
      // std::cout << "Computed transition for state " << i << "\n" << std::endl;
      progress(i,N,counter);
    }
    /* compute pre_ptr */
    abs_ptr_type sum=0;
    for(abs_type i=0; i<N; i++) {
      for(abs_type j=0; j<M; j++) {
        sum+=transition_function.m_no_pre[i*M+j];
        transition_function.m_pre_ptr[i*M+j]=sum;
      }
    }
    /* allocate memory for pre list */
    transition_function.init_transitions(T);

    /* second loop: fill pre array */
    counter=0;
    
    for(abs_type i=0; i<N; i++) {
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++) {
      /* is x an element of the overflow symbols ? */
        if(transition_function.out_of_domain[i*M+j]) 
          continue;
        /* extract lower-left and upper-bound points */
        abs_type k_lb=transition_function.corner_IDs[i*2*M+2*j];
        abs_type k_ub=transition_function.corner_IDs[i*2*M+2*j+1];
        abs_type npost=1;

        /* cell idx to coordinates */
        for(int k=dim-1; k>=0; k--) {
          /* integer coordinate of lower left corner */
          lb[k]=k_lb/NN[k];
          k_lb=k_lb-lb[k]*NN[k];
          /* integer coordinate of upper right corner */
          ub[k]=k_ub/NN[k];
          k_ub=k_ub-ub[k]*NN[k];
          /* number of grid points in each dimension in the post */
          no[k]=(ub[k]-lb[k]+1);
          /* total no of post of (i,j) */
          npost*=no[k];
          cc[k]=0;

        }

        for(abs_type k=0; k<npost; k++) {
          abs_type q=0;
          for(int l=0; l<dim; l++) 
            q+=(lb[l]+cc[l])*NN[l];
          cc[0]++;
          for(int l=0; l<dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
          /* (i,j,q) is a transition */
          transition_function.m_pre[--transition_function.m_pre_ptr[q*M+j]]=i;
        }
      }

      /* print progress */
      if(m_verbose) {
        if(counter==0)
          std::cout << "2nd loop: ";
      }
      progress(i,N,counter);
    }  
  }
 
  /** @brief get the center of cells that are used to over-approximated the
   *  attainable set associated with cell (with center x) and input u
   *
   *  @returns std::vector<abs_type>
   *  @param system_post - lambda expression as in compute_gb
   *  @param radius_post - lambda expression as in compute_gb
   *  @param x - center of cell 
   *  @param u - input
   **/
  template<class F1, class F2>
  std::vector<state_type> get_post(F1& system_post,
                                   F2& radius_post,
                                   const state_type& x,
                                 const input_type& u) const {
    /* initialize return value */
    std::vector<state_type> post {};
    /* state space dimension */
    int dim = m_state_alphabet.get_dim();
    /* variables for managing the post */
    std::vector<abs_type> lb(dim);  /* lower-left corner */
    std::vector<abs_type> ub(dim);  /* upper-right corner */
    std::vector<abs_type> no(dim);  /* number of cells per dim */
    std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
    /* get radius */
    state_type r;
    state_type eta;
    /* for out of bounds check */
    state_type lower_left;
    state_type upper_right;
    /* fill in data */
    for(int i=0; i<dim; i++) {
      eta[i]=m_state_alphabet.get_eta()[i];
      r[i]=m_state_alphabet.get_eta()[i]/2.0+m_z[i];
      lower_left[i]=m_state_alphabet.get_lower_left()[i];
      upper_right[i]=m_state_alphabet.get_upper_right()[i];
    }
    state_type xx=x;
    /* compute growth bound and numerical solution of ODE */
    radius_post(r,x,u);
    system_post(xx,u);

    /* determine post */
    abs_type npost=1;
    for(int k=0; k<dim; k++) {
      /* check for out of bounds */
      double left = xx[k]-r[k]-m_z[k];
      double right = xx[k]+r[k]+m_z[k];
      if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0) {
        return post;
      } 
      /* integer coordinate of lower left corner of post */
      lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
      /* integer coordinate of upper right corner of post */
      ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
      /* number of grid points in the post in each dimension */
      no[k]=(ub[k]-lb[k]+1);
      /* total number of post */
      npost*=no[k];
      cc[k]=0;
    }
    /* compute indices of post */
    for(abs_type k=0; k<npost; k++) {
      abs_type q=0;
      for(int l=0; l<dim; l++)  {
        q+=(lb[l]+cc[l])*m_state_alphabet.get_nn()[l];
      }
      cc[0]++;
      for(int l=0; l<dim-1; l++) {
        if(cc[l]==no[l]) {
          cc[l]=0;
          cc[l+1]++;
        }
      }
      /* (i,j,q) is a transition */    
      m_state_alphabet.itox(q,xx);
      post.push_back(xx);
    }
    return post;
  }


  /** @brief print the center of cells that are stored in the transition
   *  function as post of the cell with center x and input u
   *
   *  @param transition_function - the transition function of the abstraction
   *  @param x - center of cell 
   *  @param u - input
   **/
  void print_post(const TransitionFunction& transition_function,
                  const state_type& x,
                  const input_type& u) const {
    std::vector<state_type> post = get_post(transition_function, x, u);
    if(!post.size())
      std::cout << "\nPost is out of domain\n";
    std::cout << "\nPost states: \n";
    for(abs_type v=0; v<post.size(); v++) {
      for(int i=0; i<m_state_alphabet.get_dim(); i++) {
        std::cout << post[v][i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  /** @brief print the center of cells that are used to over-approximated the
   *  attainable set (computed according to the growth bound)
   *  associated with the cell with center x and input u
   *
   *  @param system_post - lambda expression as in compute_gb
   *  @param radius_post - lambda expression as in compute_gb
   *  @param x - center of cell 
   *  @param u - input
   **/
  template<class F1, class F2>
  void print_post_gb(F1& system_post,
                     F2& radius_post,
                     const state_type& x,
                     const input_type& u) const {
    std::vector<state_type> post = get_post(system_post,radius_post,x,u);
    if(!post.size()) {
      std::cout << "\nPost is out of domain\n";
      return;
    }

    std::cout << "\nPost states: \n";
    for(abs_type v=0; v<post.size(); v++) {
      for(int i=0; i<m_state_alphabet.get_dim(); i++) {
        std::cout << post[v][i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  /** @brief set the measurement error bound **/
  void set_measurement_error_bound(const state_type& error_bound) {
    for(int i=0; i<m_state_alphabet.get_dim(); i++) {
      m_z[i]=error_bound[i];
    }
  }
  /** @brief get measurement error bound **/
  std::vector<double> get_measruement_error_bound() {
    std::vector<double> z;
    for(int i=0; i<m_state_alphabet.get_dim(); i++) {
      z.push_back(m_z[i]);
    }
    return z;
  }
  /** @brief activate console output **/
  void verbose_on() {
    m_verbose=true;
  }
  /** @brief deactivate console output **/
  void verbose_off() {
    m_verbose=false;
  }


/*recompute transition locally*/
template<class F2, class F3, class F4=decltype(params::avoid_abs)>
  void recompute_gb(TransitionFunction& new_transition,
                    const TransitionFunction& old_transition, 
                   
                    F2& d_lb,
                    F2& d_ub,
                    F3& rs_repost,
                    F4& avoid=params::avoid_abs){

    /* number of cells */
    abs_type N=m_state_alphabet.size(); 
    /* number of inputs */
    abs_type M=m_input_alphabet.size();
    /* number of transitions (to be computed) */
    abs_ptr_type T=0; 
    /* state space dimension */
    int dim=m_state_alphabet.get_dim();
    /* for display purpose */
    abs_type counter=0;
    /* some grid information */
    std::vector<abs_type> NN=m_state_alphabet.get_nn();
    /* variables for managing the region */
    std::vector<abs_type> lb(dim);  /* lower-left corner */
    std::vector<abs_type> ub(dim);  /* upper-right corner */
    std::vector<abs_type> no(dim);  /* number of cells per dim */
    std::vector<abs_type> cc(dim);  /* coordinate of current cell in the region */
    /* radius of hyper interval containing the attainable set */
    state_type eta;
    state_type r;
    /* state and input variables */
    state_type x;
    input_type u;
    ds_type y;
    /* for out of bounds check */
    state_type lower_left;
    state_type upper_right;
    // state_type re_lb;
    // state_type re_ub;
    /* copy data from m_state_alphabet */
    for(int i=0; i<dim; i++) {
      eta[i]=m_state_alphabet.get_eta()[i];
      lower_left[i]=m_state_alphabet.get_lower_left()[i];
      upper_right[i]=m_state_alphabet.get_upper_right()[i];
      // re_lb[i]=std::max(d_lb[i]-distance[i],lower_left[i]);
      // re_ub[i]=std::min(d_ub[i]+distance[i],upper_right[i]);
    }
    /* init in transition_function the members no_pre, no_post, pre_ptr */ 
    new_transition.init_infrastructure(N,M);
    /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
    /* make if a states recompute or not.*/
    std::unique_ptr<bool[]> recomputed_mark(new bool[N*M]());
    std::unique_ptr<bool[]> input_todo(new bool[N*M]());
    /*contain the states which need to recompute*/
    std::queue<abs_type> recompute_queue; 
    std::vector<bool> out_of_region(N, true); 
    /*some variables for computing neighbours*/
    std::vector<int> tmp;
    abs_type neighbours;
    std::vector<std::vector<int>>  cc_neighbours;

    int num[]={-1,0,1};
  /*init a queue which contains the states around (d_lb, d_ub), recompute transtions for queue, 
  * if in ODE, x gets into the region(lb,ub), then enqueue the  neighbours. repeat this process untill 
  * queue become empty. make all states in queue recomputed.
  */
    /*initialize recompute_queue with the bigger box contains around states and (d_lb, d_ub)*/
    abs_type nNRegion=1;
    
   for (int i = 0; i < dim; ++i)
    {
      tmp.push_back(num[0]);
    }
    cc_neighbours.push_back(tmp);

    for (int i = dim-1; i>=0; i--)
    {
      int ncc=cc_neighbours.size();
      for (int k = 1; k < 3; ++k)
      {
        for (int j = 0; j < ncc; ++j)
        {
          tmp=cc_neighbours[j];
          tmp[i]=num[k];
          cc_neighbours.push_back(tmp);
        }
      }
    }
    int ncc = cc_neighbours.size();

    for(int k=0; k<dim; k++) {
      /* check for out of bounds */
      double left = d_lb[k]-eta[k]-m_z[k];
      double right = d_ub[k]+eta[k]+m_z[k];
      if(left <= lower_left[k]-eta[k]/2.0)
        left=lower_left[k];
      if(right >= upper_right[k]+eta[k]/2.0)
        right=upper_right[k];
      
      /* integer coordinate of lower left corner of region */
      lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
      /* integer coordinate of upper right corner of region */
      ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
      /* number of grid points in the region in each dimension */
      no[k]=(ub[k]-lb[k]+1);
      /* total number of region */
      nNRegion*=no[k];
      cc[k]=0;
    }
   
    /* compute indices of Region */
    for(abs_type k=0; k<nNRegion; k++) {
      abs_type q=0;
      for(int l=0; l<dim; l++)  {
        q+=(lb[l]+cc[l])*m_state_alphabet.get_nn()[l];
      }
      cc[0]++;
      for(int l=0; l<dim-1; l++) {
        if(cc[l]==no[l]) {
          cc[l]=0;
          cc[l+1]++;
        }
      }
      for(abs_type j=0; j<M; j++){
        if (!recomputed_mark[q*M+j]){
          recompute_queue.push(q);
          recomputed_mark[q*M+j]=true;
          input_todo[q*M+j]=true;
         
        }
      }     
      if(region(q,d_lb,d_ub,dim,eta))
        out_of_region[q]=false;  
    }
    std::cout<<"initial recomputing queue size:"<<recompute_queue.size()<<std::endl;
    abs_type coun=0;
   
   abs_type conn=0;

   /*start big loop untill the recompute_queue become empty*/
    while(!recompute_queue.empty())
    {
    /* q is a state which needs to recompute transitions, post of q*/
      abs_type q = recompute_queue.front();
      recompute_queue.pop();
      coun++;
          /* is q an element of the avoid symbols ? */
      if(avoid(q)) {
        for(abs_type j=0; j<M; j++) {
          new_transition.out_of_domain[q*M+j]=true;
        }
        continue;
      }
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++){
        //new_transition.out_of_domain[q*M+j]=false;
        if(input_todo[q*M+j]){
          input_todo[q*M+j]=false;
          /* get center x of cell */
          m_state_alphabet.itox(q,x);
          //if(x[0]<=5.92 || x[1]<=3.92 || x[0]>= 9.28 || x[1]>= 7.28 )
          //    std::cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
          /* cell radius (including measurement errors) */
         for(int k=0; k<dim; k++){
            r[k]=eta[k]/2.0+m_z[k];
            y[k]=x[k];
            y[k+dim]=r[k];
          }
    
          /* current input */
          m_input_alphabet.itox(j,u);
          /* integrate system and radius growth bound */
          /* the result is stored in x and r */
          bool intersection_with_region = false;
          rs_repost(y,u,intersection_with_region); //todo
          //if(ignore==true)
          //  continue;
          /*enqueue more neighbours of q, if q is out of new disturbance region and its trajectory has a intersection with this region*/
          for (int k = 0; k<dim; ++k)
          {
            x[k] = y[k];
            r[k] = y[k+dim];
          }
          /* determine the cells which intersect with the attainable set: 
           * discrete hyper interval of cell indices 
           * [lb[0]; ub[0]] x .... x [lb[dim-1]; ub[dim-1]]
           * covers attainable set 
           */
        
          abs_type npost=1;
          for(int k=0; k<dim; k++) {
            /* check for out of bounds */
            double left = x[k]-r[k]-m_z[k];
            double right = x[k]+r[k]+m_z[k];
            if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0){
              new_transition.out_of_domain[q*M+j]=true;
              break;
            } 

            /* integer coordinate of lower left corner of post */
            lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
            /* integer coordinate of upper right corner of post */
            ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
            /* number of grid points in the post in each dimension */
            no[k]=(ub[k]-lb[k]+1);
            /* total number of post */
            npost*=no[k];
            cc[k]=0;
          }
          //corner_IDs[q*(2*M)+2*j]=0;
          //corner_IDs[q*(2*M)+2*j+1]=0;
          if(new_transition.out_of_domain[q*M+j]){
            continue;
          }
            
          /* compute indices of post */
          for(abs_type k=0; k<npost; k++) {
            abs_type p=0;
            for(int l=0; l<dim; l++) 
              p+=(lb[l]+cc[l])*NN[l];
            cc[0]++;
            for(int l=0; l<dim-1; l++) {
              if(cc[l]==no[l]) {
                cc[l]=0;
                cc[l+1]++;
              }
            }
            /* (q,j,p) is a transition */    
            /* increment number of pres for (p,j) */ 
           // new_transition.m_no_pre[p*M+j]++;
            /* store id's of lower-left and upper-right cell */
            if(k==0)
              new_transition.corner_IDs[q*(2*M)+2*j]=p;
            if(k==npost-1)
              new_transition.corner_IDs[q*(2*M)+2*j+1]=p;
          }
            /* increment number of transitions by number of post */
           T+=npost;
           new_transition.m_no_post[q*M+j]=npost;

          if(out_of_region[q] && intersection_with_region){
            //q_neighour(q, dim, recompute_queue,recomputed_mark,NN, N);
            conn++;
            for (int i = 0; i < ncc; ++i)
            {
              neighbours=q;
              for (int k = 0; k<dim; ++k)
              {
                neighbours += cc_neighbours[i][k]*NN[k];
              }
              if(neighbours < N && !recomputed_mark[neighbours*M+j]){
                recompute_queue.push(neighbours);
                recomputed_mark[neighbours*M+j]=true;
                input_todo[neighbours*M+j]=true;        
             
              }
            }
          }  
        }
      }
    }
  
    std::cout<<"recomputing transitions untill queue empty";
    std::cout<<"total number of recomputing states * inputs:" <<coun<<std::endl;
      /*copy from old transtions*/
      for(abs_type i = 0; i < N; ++i )
      { 
        for (abs_type j = 0; j < M; ++j)
        {
          // if(recomputed_mark[i*M+j]){
          //   if(new_transition.out_of_domain[i*M+j]!=standard_transition.out_of_domain[i*M+j])
          //     std::cout<<"error"<<i<<std::endl;
          // }
          //new_transition.out_of_domain[i*M+j]=standard_transition.out_of_domain[i*M+j];
          if(!recomputed_mark[i*M+j]){
            new_transition.out_of_domain[i*M+j]=old_transition.out_of_domain[i*M+j];
            new_transition.m_no_post[i*M+j]=old_transition.m_no_post[i*M+j];
            T+=new_transition.m_no_post[i*M+j];
            new_transition.corner_IDs[i*(2*M)+2*j]=old_transition.corner_IDs[i*(2*M)+2*j];
            new_transition.corner_IDs[i*(2*M)+2*j+1]=old_transition.corner_IDs[i*(2*M)+2*j+1];
          }           
        }
      }

       counter=0;
       int n_lb=0;
     for(abs_type i=0; i<N; i++) {
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++) {
      /* is x an element of the overflow symbols ? */
        if(new_transition.out_of_domain[i*M+j]) 
          continue;
        /* extract lower-left and upper-bound points */
        abs_type k_lb=new_transition.corner_IDs[i*(2*M)+2*j];
        if(k_lb==0)
          n_lb++;
        abs_type k_ub=new_transition.corner_IDs[i*2*M+2*j+1];
        abs_type npost=1;

        /* cell idx to coordinates */
        for(int k=dim-1; k>=0; k--) {
          /* integer coordinate of lower left corner */
          lb[k]=k_lb/NN[k];
          k_lb=k_lb-lb[k]*NN[k];
          /* integer coordinate of upper right corner */
          ub[k]=k_ub/NN[k];
          k_ub=k_ub-ub[k]*NN[k];
          /* number of grid points in each dimension in the post */
          no[k]=(ub[k]-lb[k]+1);
          /* total no of post of (i,j) */
          npost*=no[k];
          cc[k]=0;
        }

        for(abs_type k=0; k<npost; k++) {
          abs_type q=0;
          for(int l=0; l<dim; l++) 
            q+=(lb[l]+cc[l])*NN[l];
          cc[0]++;
          for(int l=0; l<dim-1; l++) {
            if(cc[l]==no[l]) {
              cc[l]=0;
              cc[l+1]++;
            }
          }
          /* (i,j,q) is a transition */
          new_transition.m_no_pre[q*M+j]++;
        }
      }
        /* print progress */
        if(m_verbose) {
          if(counter==0)
            std::cout << "2nd loop: ";
        }
        progress(i,N,counter);
      }

      /* compute pre_ptr */
      abs_ptr_type sum=0;
      for(abs_type i=0; i<N; i++) {
        for(abs_type j=0; j<M; j++) {
          sum+=new_transition.m_no_pre[i*M+j];
          new_transition.m_pre_ptr[i*M+j]=sum; 
          // if (new_transition.out_of_domain[i*M+j]!=standard_transition.out_of_domain[i*M+j])
          //     {
          //       std::cout<<i<<" "<<j<<recomputed_mark[i*M+j]<<" "<<new_transition.out_of_domain[i*M+j]<<" "<<standard_transition.out_of_domain[i*M+j]<<std::endl;
          //     }    
        }
      }
        
      /* allocate memory for pre list */
      new_transition.init_transitions(T);
      
      /*loop: fill pre array */
      counter=0;
      for(abs_type i=0; i<N; i++) {
        /* loop over all inputs */
        for(abs_type j=0; j<M; j++) {
        /* is x an element of the overflow symbols ? */
          if(new_transition.out_of_domain[i*M+j]) 
            continue;
          /* extract lower-left and upper-bound points */
          abs_type k_lb=new_transition.corner_IDs[i*(2*M)+2*j];
          abs_type k_ub=new_transition.corner_IDs[i*(2*M)+2*j+1];
          abs_type npost=1;

          /* cell idx to coordinates */
          for(int k=dim-1; k>=0; k--) {
            /* integer coordinate of lower left corner */
            lb[k]=k_lb/NN[k];
            k_lb=k_lb-lb[k]*NN[k];
            /* integer coordinate of upper right corner */
            ub[k]=k_ub/NN[k];
            k_ub=k_ub-ub[k]*NN[k];
            /* number of grid points in each dimension in the post */
            no[k]=(ub[k]-lb[k]+1);
            /* total no of post of (i,j) */
            npost*=no[k];
            cc[k]=0;
          }

          for(abs_type k=0; k<npost; k++) {
            abs_type p=0;
            for(int l=0; l<dim; l++) 
              p+=(lb[l]+cc[l])*NN[l];
            cc[0]++;
            for(int l=0; l<dim-1; l++) {
              if(cc[l]==no[l]) {
                cc[l]=0;
                cc[l+1]++;
              }
            }
            /* (i,j,p) is a transition */
              new_transition.m_pre[--new_transition.m_pre_ptr[p*M+j]]=i;
          }
        }
        /* print progress */
        if(m_verbose) {
          if(counter==0)
            std::cout << "3nd loop: ";
        }
        progress(i,N,counter);
      }
   
    write_to_filebb(m_state_alphabet,m_input_alphabet,recomputed_mark,"recomputation_lazy");
  /*for all states, if recomputed_mark==ture, to get the differences of old. 
  * if the recomputed_mark==false, to copy the transitions from old transitions*/

    // for (abs_type i = 0; i < N; ++i)
    // {
    //   if (recomputed_mark[i]==true)
    //     std::queue<abs_type> dif = get_difference(i);//todo
    //     diff_queue.push_back(dif); //tddo
    // }
  
  }//function closed
 
bool region(abs_type q,state_type d_lb,state_type d_ub, int dim, state_type eta){
  state_type x;
  m_state_alphabet.itox(q,x);
  bool belong=true;
  for(int i=0; i<dim; i++){
    if (d_lb[i]-eta[i]/2.0 > x[i] || x[i] > d_ub[i]+eta[i]/2.0)
      belong = false;
  }
  return belong;
}

/*recompute transition locally*/
template<class F1, class F2, class F3, class F4=decltype(params::avoid_abs)>
  void recompute_mr(TransitionFunction& new_transition_com,
                    const TransitionFunction& old_transition,
                    F1& distance,
                    F2& d_lb,
                    F2& d_ub,
                    F3& rs_post,
                    F4& avoid=params::avoid_abs){

    /* number of cells */
    abs_type N=m_state_alphabet.size(); 
    /* number of inputs */
    abs_type M=m_input_alphabet.size();
    /* number of transitions (to be computed) */
    abs_ptr_type T=0; 
    /* state space dimension */
    int dim=m_state_alphabet.get_dim();
    /* for display purpose */
    abs_type counter=0;
    /* some grid information */
    std::vector<abs_type> NN=m_state_alphabet.get_nn();
    /* variables for managing the region */
    std::vector<abs_type> lb(dim);  /* lower-left corner */
    std::vector<abs_type> ub(dim);  /* upper-right corner */
    std::vector<abs_type> no(dim);  /* number of cells per dim */
    std::vector<abs_type> cc(dim);  /* coordinate of current cell in the region */
    /* radius of hyper interval containing the attainable set */
    state_type eta;
    state_type r;
    /* state and input variables */
    state_type x;
    input_type u;
    ds_type y;
    /* for out of bounds check */
    state_type lower_left;
    state_type upper_right;
    state_type re_lb;
    state_type re_ub;
    /* copy data from m_state_alphabet */
    for(int i=0; i<dim; i++) {
      eta[i]=m_state_alphabet.get_eta()[i];
      lower_left[i]=m_state_alphabet.get_lower_left()[i];
      upper_right[i]=m_state_alphabet.get_upper_right()[i];
      re_lb[i]=std::max(d_lb[i]-distance[i],lower_left[i]);
      re_ub[i]=std::min(d_ub[i]+distance[i],upper_right[i]);
    }
    /* init in transition_function the members no_pre, no_post, pre_ptr */ 
    new_transition_com.init_infrastructure(N,M);
    /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
    /* make if a states recompute or not.*/
    std::vector<bool> recomputed_mark(N, false);
    std::unique_ptr<bool[]> input_todo(new bool[N*M]());
    /*contain the states which need to recompute*/
    std::queue<abs_type> recompute_queue; 
    std::vector<bool> out_of_region(N, true); 

  /*init a queue which contains the states around (d_lb, d_ub), recompute transtions for queue, 
  * if in ODE, x gets into the region(lb,ub), then enqueue the  neighbours. repeat this process untill 
  * queue become empty. make all states in queue recomputed.
  */
    /*initialize recompute_queue with the bigger box contains around states and (d_lb, d_ub)*/
    abs_type nNRegion=1;
    

    for(int k=0; k<dim; k++) {
      /* check for out of bounds */
      double left = re_lb[k]-m_z[k];
      double right = re_ub[k]+m_z[k];
      if(left <= lower_left[k]-eta[k]/2.0)
        left=lower_left[k];
      if(right >= upper_right[k]+eta[k]/2.0)
        right=upper_right[k];
      
      /* integer coordinate of lower left corner of region */
      lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
      /* integer coordinate of upper right corner of region */
      ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
      /* number of grid points in the region in each dimension */
      no[k]=(ub[k]-lb[k]+1);
      /* total number of region */
      nNRegion*=no[k];
      cc[k]=0;
    }
   
    /* compute indices of Region */
    for(abs_type k=0; k<nNRegion; k++) {
      abs_type q=0;
      for(int l=0; l<dim; l++)  {
        q+=(lb[l]+cc[l])*m_state_alphabet.get_nn()[l];
      }
      cc[0]++;
      for(int l=0; l<dim-1; l++) {
        if(cc[l]==no[l]) {
          cc[l]=0;
          cc[l+1]++;
        }
      }

      if (!recomputed_mark[q]){
        recompute_queue.push(q);
        recomputed_mark[q]=true;
      }
    }
    std::cout<<"initial recomputing queue size:"<<recompute_queue.size()<<std::endl;
    abs_type coun=0;

   
   abs_type conn=0;
  
   /*start big loop untill the recompute_queue become empty*/
    while(!recompute_queue.empty())
    {
    /* q is a state which needs to recompute transitions, post of q*/
      abs_type q = recompute_queue.front();
      recompute_queue.pop();
      coun++;
          /* is q an element of the avoid symbols ? */
      if(avoid(q)) {
        for(abs_type j=0; j<M; j++) {
          new_transition_com.out_of_domain[q*M+j]=true;
        }
        continue;
      }
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++) {
        new_transition_com.out_of_domain[q*M+j]=false;         
          /* get center x of cell */
          m_state_alphabet.itox(q,x);
          //if(x[0]<=5.92 || x[1]<=3.92 || x[0]>= 9.28 || x[1]>= 7.28 )
          //    std::cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
          /* cell radius (including measurement errors) */
         for(int k=0; k<dim; k++){
            r[k]=eta[k]/2.0+m_z[k];
            y[k]=x[k];
            y[k+dim]=r[k];
          }
        
          /* current input */
          m_input_alphabet.itox(j,u);
          /* integrate system and radius growth bound */
          /* the result is stored in x and r */

          rs_post(y,u); //todo
          
          /*enqueue more neighbours of q, if q is out of new disturbance region and its trajectory has a intersection with this region*/
          for (int k = 0; k<dim; ++k)
          {
            x[k] = y[k];
            r[k] = y[k+dim];
          }
          /* determine the cells which intersect with the attainable set: 
           * discrete hyper interval of cell indices 
           * [lb[0]; ub[0]] x .... x [lb[dim-1]; ub[dim-1]]
           * covers attainable set 
           */
         
          abs_type npost=1;
          for(int k=0; k<dim; k++) {
            /* check for out of bounds */
            double left = x[k]-r[k]-m_z[k];
            double right = x[k]+r[k]+m_z[k];
            if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0)  {
              new_transition_com.out_of_domain[q*M+j]=true;
              break;
            } 

            /* integer coordinate of lower left corner of post */
            lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
            /* integer coordinate of upper right corner of post */
            ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
            /* number of grid points in the post in each dimension */
            no[k]=(ub[k]-lb[k]+1);
            /* total number of post */
            npost*=no[k];
            cc[k]=0;
          }
          //corner_IDs[q*(2*M)+2*j]=0;
          //corner_IDs[q*(2*M)+2*j+1]=0;
          if(new_transition_com.out_of_domain[q*M+j])
            continue;
            
          /* compute indices of post */
          for(abs_type k=0; k<npost; k++) {
            abs_type p=0;
            for(int l=0; l<dim; l++) 
              p+=(lb[l]+cc[l])*NN[l];
            cc[0]++;
            for(int l=0; l<dim-1; l++) {
              if(cc[l]==no[l]) {
                cc[l]=0;
                cc[l+1]++;
              }
            }
            /* (q,j,p) is a transition */    
            /* increment number of pres for (p,j) */ 
           // new_transition.m_no_pre[p*M+j]++;
            /* store id's of lower-left and upper-right cell */
            if(k==0)
              new_transition_com.corner_IDs[q*(2*M)+2*j]=p;
            if(k==npost-1)
              new_transition_com.corner_IDs[q*(2*M)+2*j+1]=p;
            }
            /* increment number of transitions by number of post */
           T+=npost;
           new_transition_com.m_no_post[q*M+j]=npost;               
      }
    }
    std::cout<<"recomputing transitions untill queue empty";
  
      
      std::cout<<"total number of recomputing states * inputs:" <<coun<<"\n"<<conn<<std::endl;
      /*copy from old transtions*/
      for (abs_type i = 0; i < N; ++i)
      { 
        for (abs_type j = 0; j < M; ++j)
        {
          if(!recomputed_mark[i]){
            new_transition_com.out_of_domain[i*M+j]=old_transition.out_of_domain[i*M+j];
            // new_transition.m_no_pre[i*M+j]=old_transition.m_no_pre[i*M+j];
            new_transition_com.m_no_post[i*M+j]=old_transition.m_no_post[i*M+j];
            T+=new_transition_com.m_no_post[i*M+j];
            new_transition_com.corner_IDs[i*(2*M)+2*j]=old_transition.corner_IDs[i*(2*M)+2*j];
            new_transition_com.corner_IDs[i*(2*M)+2*j+1]=old_transition.corner_IDs[i*(2*M)+2*j+1];
          }

          // if(new_transition.corner_IDs[i*(2*M)+2*j]!=standard_transition.corner_IDs[i*(2*M)+2*j] ||
          //   new_transition.corner_IDs[i*(2*M)+2*j+1]!=standard_transition.corner_IDs[i*(2*M)+2*j+1])
          //   std::cout<<"here:"<<i<<" "<<j<<" "<<recomputed_mark[i*M+j]<<std::endl;
           
        }
      }

       counter=0;
      for(abs_type i=0; i<N; i++) {
        /* loop over all inputs */
        for(abs_type j=0; j<M; j++) {
        /* is x an element of the overflow symbols ? */
          if(new_transition_com.out_of_domain[i*M+j]) 
            continue;
          /* extract lower-left and upper-bound points */
          abs_type k_lb=new_transition_com.corner_IDs[i*2*M+2*j];
          abs_type k_ub=new_transition_com.corner_IDs[i*2*M+2*j+1];
          abs_type npost=1;

          /* cell idx to coordinates */
          for(int k=dim-1; k>=0; k--) {
            /* integer coordinate of lower left corner */
            lb[k]=k_lb/NN[k];
            k_lb=k_lb-lb[k]*NN[k];
            /* integer coordinate of upper right corner */
            ub[k]=k_ub/NN[k];
            k_ub=k_ub-ub[k]*NN[k];
            /* number of grid points in each dimension in the post */
            no[k]=(ub[k]-lb[k]+1);
            /* total no of post of (i,j) */
            npost*=no[k];
            cc[k]=0;

          }

          for(abs_type k=0; k<npost; k++) {
            abs_type p=0;
            for(int l=0; l<dim; l++) 
              p+=(lb[l]+cc[l])*NN[l];
            cc[0]++;
            for(int l=0; l<dim-1; l++) {
              if(cc[l]==no[l]) {
                cc[l]=0;
                cc[l+1]++;
              }
            }
            /* (i,j,p) is a transition */
            new_transition_com.m_no_pre[p*M+j]++;
          }
        }
        /* print progress */
        if(m_verbose) {
          if(counter==0)
            std::cout << "2nd loop: ";
        }
        progress(i,N,counter);
      }

      /* compute pre_ptr */
      abs_ptr_type sum=0;
      for(abs_type i=0; i<N; i++) {
        for(abs_type j=0; j<M; j++) {
          sum+=new_transition_com.m_no_pre[i*M+j];
          new_transition_com.m_pre_ptr[i*M+j]=sum;
        }
      }
      
      
      /* allocate memory for pre list */
      new_transition_com.init_transitions(T);
      
      /*loop: fill pre array */
      counter=0;
      for(abs_type i=0; i<N; i++) {
        /* loop over all inputs */
        for(abs_type j=0; j<M; j++) {
        /* is x an element of the overflow symbols ? */
          if(new_transition_com.out_of_domain[i*M+j]) 
            continue;
          /* extract lower-left and upper-bound points */
          abs_type k_lb=new_transition_com.corner_IDs[i*2*M+2*j];
          abs_type k_ub=new_transition_com.corner_IDs[i*2*M+2*j+1];
          abs_type npost=1;

          /* cell idx to coordinates */
          for(int k=dim-1; k>=0; k--) {
            /* integer coordinate of lower left corner */
            lb[k]=k_lb/NN[k];
            k_lb=k_lb-lb[k]*NN[k];
            /* integer coordinate of upper right corner */
            ub[k]=k_ub/NN[k];
            k_ub=k_ub-ub[k]*NN[k];
            /* number of grid points in each dimension in the post */
            no[k]=(ub[k]-lb[k]+1);
            /* total no of post of (i,j) */
            npost*=no[k];
            cc[k]=0;
          }

          for(abs_type k=0; k<npost; k++) {
            abs_type p=0;
            for(int l=0; l<dim; l++) 
              p+=(lb[l]+cc[l])*NN[l];
            cc[0]++;
            for(int l=0; l<dim-1; l++) {
              if(cc[l]==no[l]) {
                cc[l]=0;
                cc[l+1]++;
              }
            }
            /* (i,j,p) is a transition */
            new_transition_com.m_pre[--new_transition_com.m_pre_ptr[p*M+j]]=i;
          }
        }
        /* print progress */
        if(m_verbose) {
          if(counter==0)
            std::cout << "3nd loop: ";
        }
        progress(i,N,counter);
      }
      write_to_fileb(m_state_alphabet,recomputed_mark,"recomputation_max");

    } 
};

} /* close namespace */
#endif /* ABSTRACTION_HH_ */
