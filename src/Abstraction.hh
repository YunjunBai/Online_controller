/*
 * Abstraction.hh
 *
 *  
 */

/** @file **/
#ifndef ABSTRACTION_HH_
#define ABSTRACTION_HH_

#include <iostream>
#include <cstring>
#include <memory>
#include <vector>
#include <array>
#include <tuple>
#include <chrono>
#include <thread>
#include <atomic>
#include "UniformGrid.hh"
#include "TransitionFunction.hh"
#include "Disturbance.hh"
#include "InputOutput.hh"
#include "TicToc.hh"

#define CHUNK_SIZE 800

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
    template<std::size_t state_dim, std::size_t input_dim>
        class Abstraction:public TransitionFunction {

            using state_type = std::array<double,state_dim>;
            using disturbance_type = std::array<double, state_dim>;
            using input_type = std::array<double,input_dim>;
            using ds_type = std::array<double, 2*state_dim>;

            private:
                /* grid information of state alphabet */
                const UniformGrid m_state_alphabet;
                /* grid information of input alphabet */
                const UniformGrid m_input_alphabet;
                /* measurement error bound */
                std::array<double,state_dim> m_z;
                /* print progress to the console (default m_verbose=true) */
                bool m_verbose=true;

                /* number of cells */
                const abs_type N; 
                /* number of inputs */
                const abs_type M; 
                /* some grid information */
                const std::vector<abs_type> NN;
                /* radius of hyper interval containing the attainable set */
                state_type eta;
                /* for out of bounds check */
                state_type lower_left;
                state_type upper_right;

                /* to display the progress of the computation of the abstraction */
                void progress(const abs_type& i, abs_type& counter) {
                    if(!m_verbose)
                        return;
                    if(((double)i/(double)N*100)>=counter){
                        if((counter%10)==0)
                            std::cout << counter ;
                        else if((counter%2)==0) 
                            std::cout << ".";
                        counter++;
                        //std::cout << "thread(" << omp_get_thread_num() << "): " ;
                    }
                    std::flush(std::cout); 
                    if(i==(N-1))
                        std::cout <<"100\n";

                }
                void progress_re(abs_type& counter, abs_type& counter_states) {
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

                void fill_pre(TransitionFunction& transition_function)
                {
                    /* for display purpose */
                    abs_type counter=0;
                    /* compute pre_ptr */
                    abs_ptr_type sum=0;
                    for(abs_type i=0; i<N; i++) {
                        for(abs_type j=0; j<M; j++) {
                            sum+=transition_function.m_no_pre[i*M+j];
                            transition_function.m_pre_ptr[i*M+j]=sum;
                        }
                    }
                    std::cout<<"sum:"<<sum<<std::endl;
                    /* allocate memory for pre list */
                    transition_function.init_transitions(sum);
#pragma omp parallel
                    {
                        /* variables for managing the region */
                        std::array<abs_type,state_dim> lb;  /* lower-left corner */
                        std::array<abs_type,state_dim> ub;  /* upper-right corner */
                        std::array<abs_type,state_dim> no;  /* number of cells per dim */
                        std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */
                        /* second loop: fill pre array */
#pragma omp for schedule(dynamic, CHUNK_SIZE)
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
                                for(int k=state_dim-1; k>=0; k--) {
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
                                    for(std::size_t l=0; l<state_dim; l++) 
                                        q+=(lb[l]+cc[l])*NN[l];
                                    cc[0]++;
                                    for(std::size_t l=0; l<state_dim-1; l++) {
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
                                    std::cout << "fill pre loop: ";
                            }
                            progress(i,counter);
                        } // end for
                    } // end parallel
                }

                inline void compute_post_indices(
                    TransitionFunction& transition_function,
                    const abs_type i,
                    const abs_type j,
                    const abs_type npost,

                    std::array<abs_type,state_dim>  lb,
                    std::array<abs_type,state_dim> no,
                    std::array<abs_type,state_dim> & cc)

                {
                   // std::cout<<"excuse this function"<<std::endl;
                    /* compute indices of post */
                    for(abs_type k=0; k<npost; k++) {
                        abs_type q=0;
                        for(std::size_t l=0; l<state_dim; l++) 
                            q+=(lb[l]+cc[l])*NN[l];
                        cc[0]++;
                        for(std::size_t l=0; l<state_dim-1; l++) {
                            if(cc[l]==no[l]) {
                                cc[l]=0;
                                cc[l+1]++;
                            }
                        }
                        /* (i,j,q) is a transition */    
                        /* increment number of pres for (q,j) */ 

                        abs_type index = q*M+j;
                        #pragma omp atomic
                        transition_function.m_no_pre[index]++;

                        /* store id's of lower-left and upper-right cell */
                        if(k==0)
                            transition_function.corner_IDs[i*(2*M)+2*j]=q;
                        if(k==npost-1)
                            transition_function.corner_IDs[i*(2*M)+2*j+1]=q;
                    }
                }
                                    
                inline void process_post(
                    TransitionFunction & transition_function,
                    const abs_type i,
                    const abs_type j,
                    state_type & x,
                    ds_type & y,
                    state_type & r,
                    abs_type & npost,
                    std::array<abs_type,state_dim> & lb,
                    std::array<abs_type,state_dim> & ub,
                    std::array<abs_type,state_dim> & no,
                    std::array<abs_type,state_dim> & cc)
                {
                    for (std::size_t k = 0; k<state_dim; ++k)
                    {
                        x[k]=y[k];
                        r[k]=y[k+state_dim];
                    }
                    /* determine the cells which intersect with the attainable set: 
                     * discrete hyper interval of cell indices 
                     * [lb[0]; ub[0]] x .... x [lb[state_dim-1]; ub[dim-1]]
                     * covers attainable set 
                     */
                    for(std::size_t k=0; k<state_dim; k++) {
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
                    m_z(),
                    N(state_alphabet.size()),
                    M(input_alphabet.size()),
                    NN(state_alphabet.get_nn())
                    {
                        /* default value of the measurement error 
                         * (heurisitc to prevent rounding errors)*/
                        for(std::size_t i=0; i<state_dim; i++)
                            m_z[i]=m_state_alphabet.get_eta()[i]/1e10;
                        /* copy data from m_state_alphabet */
                        for(std::size_t i=0; i<state_dim; i++) {
                            eta[i]=m_state_alphabet.get_eta()[i];
                            lower_left[i]=m_state_alphabet.get_lower_left()[i];
                            upper_right[i]=m_state_alphabet.get_upper_right()[i];
                        }
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

                        /* for display purpose */
                        abs_type counter=0;
                        /* init in transition_function the members no_pre, no_post, pre_ptr */ 
                        transition_function.init_infrastructure(N,M);
                        /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */

                        /*
                         * first loop: compute corner_IDs:
                         * corner_IDs[i*M+j][0] = lower-left cell index of over-approximation of attainable set 
                         * corner_IDs[i*M+j][1] = upper-right cell index of over-approximation of attainable set 
                         */
                        /* loop over all cells */
#pragma omp parallel shared(counter)
                        {
                            //std::cout<<"the number of threads:"<<omp_get_num_threads()<<std::endl;
                            state_type x;
                            input_type u;
                            ds_type y;
                            state_type r; 
                            /* variables for managing the region */
                            std::array<abs_type,state_dim> lb;  /* lower-left corner */
                            std::array<abs_type,state_dim> ub;  /* upper-right corner */
                            std::array<abs_type,state_dim> no;  /* number of cells per dim */
                            std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */

#pragma omp for schedule(dynamic, CHUNK_SIZE)
                            for(abs_type i=0; i<N; i++) {

                                if(avoid(i)) {
                                    for(abs_type j=0; j<M; j++) {
                                        transition_function.out_of_domain[i*M+j]=true;
                                    }
                                    continue;
                                }
                                /* loop over all inputs */
                                for(abs_type j=0; j<M; j++) {
                                    /* is i an element of the avoid symbols ? */

                                    transition_function.out_of_domain[i*M+j]=false;
                                    /* get center x of cell */
                                    m_state_alphabet.itox(i,x);
                                    /* cell radius (including measurement errors) */
                                    for(std::size_t k=0; k<state_dim; k++){
                                        r[k]=eta[k]/2.0+m_z[k];
                                        y[k]=x[k];
                                        y[k+state_dim]=r[k];
                                    }

                                    /* current input */
                                    m_input_alphabet.itox(j,u);
                                    /* integrate system and radius growth bound */
                                    /* the result is stored in x and r */

                                    rs_post(y,u);
                                    abs_type npost=1;
                                    process_post(transition_function, i, j, x, y, r, npost, lb, ub, no, cc);
                                    //std::cout<<"nuem"<<npost<<std::endl;
                                    if(transition_function.out_of_domain[i*M+j]) 
                                        continue;

                                    compute_post_indices(transition_function, i, j, npost, lb, no, cc);
                                    /* increment number of transitions by number of post */
                                    transition_function.m_no_post[i*M+j]=npost;
                                }
                                /* print progress */
                                if(m_verbose) {
                                    if(counter==0)
                                        std::cout << "1st loop: ";
                                }
                                // std::cout << "Computed transition for state " << i << "\n" << std::endl;
                                progress(i,counter);
                            } // end for
                        } // end parallel

                        fill_pre(transition_function);
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
                        /* variables for managing the post */
                        std::array<abs_type,state_dim> lb;  /* lower-left corner */
                        std::array<abs_type,state_dim> ub;  /* upper-right corner */
                        std::array<abs_type,state_dim> no;  /* number of cells per dim */
                        std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */
                        /* get radius */
                        state_type r;
                        /* fill in data */
                        for(int i=0; i<state_dim; i++) {
                            r[i]=m_state_alphabet.get_eta()[i]/2.0+m_z[i];
                        }
                        state_type xx=x;
                        /* compute growth bound and numerical solution of ODE */
                        radius_post(r,x,u);
                        system_post(xx,u);

                        /* determine post */
                        abs_type npost=1;
                        for(int k=0; k<state_dim; k++) {
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
                            for(int l=0; l<state_dim; l++)  {
                                q+=(lb[l]+cc[l])*NN[l];
                            }
                            cc[0]++;
                            for(int l=0; l<state_dim-1; l++) {
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
                std::vector<double> get_measurement_error_bound() {
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
                            std::queue<abs_type>& diff,
                            const TransitionFunction& old_transition, 
                            F2& d_lb,
                            F2& d_ub,
                            F3& rs_repost,
                            F4& avoid=params::avoid_abs){

                        /* for display purpose */
                        abs_type counter=0;
                        /* variables for managing the region */
                        std::array<abs_type,state_dim> lb;  /* lower-left corner */
                        std::array<abs_type,state_dim> ub;  /* upper-right corner */
                        std::array<abs_type,state_dim> no;  /* number of cells per dim */
                        std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */

                        /* init in transition_function the members no_pre, no_post, pre_ptr */ 
                        new_transition.init_infrastructure(N,M);
                        /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
                        /* make if a states recompute or not.*/
                        std::unique_ptr<std::atomic_flag[]> recomputed_mark(new std::atomic_flag[N*M]());
                        for (abs_type i = 0; i < N*M; i++) {
                            recomputed_mark[i].clear(); //ATOMIC_FLAG_INIT;
                        }
                        /*contain the states which need to recompute*/
                        std::queue<std::pair<abs_type,abs_type>> recompute_queue; 
                        std::vector<bool> out_of_region(N, true); 
                        std::vector<bool> diff_done(N,false);

                        /*some variables for computing neighbours*/
                        std::vector<int> tmp;
                        std::vector<std::vector<int>>  cc_neighbours;

                        int num[]={-1,0,1};
                        /*init a queue which contains the states around (d_lb, d_ub), recompute transtions for queue, 
                         * if in ODE, x gets into the region(lb,ub), then enqueue the  neighbours. repeat this process untill 
                         * queue become empty. make all states in queue recomputed.
                         */
                        /*initialize recompute_queue with the bigger box contains around states and (d_lb, d_ub)*/
                        abs_type nNRegion=1;

                        for (std::size_t i = 0; i < state_dim; ++i)
                        {
                            tmp.push_back(num[0]);
                        }
                        cc_neighbours.push_back(tmp);

                        for (int i = state_dim-1; i>=0; i--)
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

                        for(std::size_t k=0; k<state_dim; k++) {
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
                            for(std::size_t l=0; l<state_dim; l++)  {
                                q+=(lb[l]+cc[l])*m_state_alphabet.get_nn()[l];
                            }
                            cc[0]++;
                            for(std::size_t l=0; l<state_dim-1; l++) {
                                if(cc[l]==no[l]) {
                                    cc[l]=0;
                                    cc[l+1]++;
                                }
                            }
                            for(abs_type j=0; j<M; j++){
                                bool already_done = recomputed_mark[q*M+j].test_and_set();
                                if (!already_done){
                                    recompute_queue.emplace(q,j);
                                }
                            }     
                            if(region(q,d_lb,d_ub,eta))
                                out_of_region[q]=false;  
                        }
                        std::cout<<"initial recomputing queue size:"<<recompute_queue.size()<<std::endl;


                        abs_type conn=0;
                        abs_type coun=0;
                        /*start big loop untill the recompute_queue become empty*/

                        //numbers of thread currently working
                        int workers = 0;

#pragma omp parallel
                        {
                            //std::cout<<"the number of threads:"<<omp_get_num_threads()<<std::endl;
                            state_type x;
                            input_type u;
                            ds_type y;
                            state_type r; 
                            /* variables for managing the region */
                            std::array<abs_type,state_dim> lb;  /* lower-left corner */
                            std::array<abs_type,state_dim> ub;  /* upper-right corner */
                            std::array<abs_type,state_dim> no;  /* number of cells per dim */
                            std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */
                            int chunk = 0;
                            std::array<abs_type,CHUNK_SIZE> work_list_states;
                            std::array<abs_type,CHUNK_SIZE> work_list_inputs;
                            bool working = true;
                            int new_chunk = 0;
                            std::array<abs_type,CHUNK_SIZE> new_work_list_states;
                            std::array<abs_type,CHUNK_SIZE> new_work_list_inputs;

#pragma omp atomic
                            workers++;

                            //as long as we have busy workers
                            while(workers > 0) {
                                //get a chuck of states to work on
#pragma omp critical
                                {
                                    while(chunk < CHUNK_SIZE && !recompute_queue.empty()) {
                                        work_list_states[chunk] = std::get<0>(recompute_queue.front());
                                        work_list_inputs[chunk] = std::get<1>(recompute_queue.front()); 
                                        recompute_queue.pop();
                                        chunk++;
                                    }
                                    coun += chunk;
                                }

                                if (chunk == 0) {
                                    //no more work, says you are done
                                    working = false;
#pragma omp atomic
                                    workers--;
                                    //wait a bit, someone else might still be busy and will enqueue new states
                                    std::this_thread::sleep_for(std::chrono::milliseconds(200));
                                } else {
                                    if (!working) {
                                        //back to work!
                                        working = true;
#pragma omp atomic
                                        workers++;
                                    }
                                    for (int i = 0; i < chunk; i++) {
                                        abs_type q = work_list_states[i];
                                        abs_type j = work_list_inputs[i];
                                        /* is q an element of the avoid symbols ? */
                                        if(avoid(q)) {
                                            for(abs_type j=0; j<M; j++) {
                                                new_transition.out_of_domain[q*M+j]=true;
                                            }
                                            continue;
                                        }
                                        /* get center x of cell */
                                        m_state_alphabet.itox(q,x);
                                        /* cell radius (including measurement errors) */
                                        for(std::size_t k=0; k<state_dim; k++){
                                            r[k]=eta[k]/2.0+m_z[k];
                                            y[k]=x[k];
                                            y[k+state_dim]=r[k];
                                        }

                                        /* current input */
                                        m_input_alphabet.itox(j,u);
                                        /* integrate system and radius growth bound */
                                        /* the result is stored in x and r */
                                       bool intersection_with_region = false;
                                        rs_repost(y,u,intersection_with_region); 
                                        abs_type npost=1;
                                        process_post(new_transition, q, j, x, y, r, npost, lb, ub, no, cc);

                                        if(new_transition.out_of_domain[q*M+j]){
                                            continue;
                                        }

                                        compute_post_indices(new_transition, q, j, npost, lb, no, cc);
                                        /* increment number of transitions by number of post */
                                        new_transition.m_no_post[q*M+j]=npost;

                                        /*enqueue more neighbours of q, if q is out of new disturbance region and its trajectory has a intersection with this region*/
                                        if(out_of_region[q] && intersection_with_region){
                                            conn++;
                                            for (int i = 0; i < ncc; ++i)
                                            {
                                                abs_type neighbours=q;
                                                for (std::size_t k = 0; k<state_dim; ++k)
                                                {
                                                    neighbours += cc_neighbours[i][k]*NN[k];
                                                }
                                                if(neighbours < N) {
                                                    bool already_done = recomputed_mark[neighbours*M+j].test_and_set();
                                                    if (!already_done){
                                                        new_work_list_states[new_chunk] = neighbours;
                                                        new_work_list_inputs[new_chunk] = j;
                                                        new_chunk++;
                                                        //about to overflow the new work list
                                                        if (new_chunk == CHUNK_SIZE) {
                                                            #pragma omp critical
                                                            for (int i = 0; i < new_chunk; i++) {
                                                                recompute_queue.emplace(new_work_list_states[i],new_work_list_inputs[i]);
                                                            }
                                                            new_chunk = 0;
                                                        }
                                                    }
                                                }
                                            }
                                        }  
                                    }
                                    chunk = 0;
                                  // put new stuff in the work queue
                                    #pragma omp critical
                                    for (int i = 0; i < new_chunk; i++) {
                                        recompute_queue.emplace(new_work_list_states[i],new_work_list_inputs[i]);
                                    }
                                    new_chunk = 0;
                                }
                            }
                        } // end parallel
                        
                        std::unique_ptr<bool[]> recomputed(new bool[N*M]());
                        for (abs_type i = 0; i < N*M; i++) {
                            recomputed[i] = recomputed_mark[i].test_and_set();
                        }
                        std::cout<<"total " <<coun<<std::endl;
                        /*copy from old transtions*/
                        for(abs_type i = 0; i < N; ++i )
                        { 
                            for (abs_type j = 0; j < M; ++j)
                            {
                                if(recomputed[i*M+j] && !diff_done[i] && !avoid(i)){
                                    diff.push(i);
                                    diff_done[i]=true;
                                }
                                //new_transition.out_of_domain[i*M+j]=standard_transition.out_of_domain[i*M+j];
                                else{
                                    new_transition.out_of_domain[i*M+j]=old_transition.out_of_domain[i*M+j];
                                    new_transition.m_no_post[i*M+j]=old_transition.m_no_post[i*M+j];
                                    new_transition.corner_IDs[i*(2*M)+2*j]=old_transition.corner_IDs[i*(2*M)+2*j];
                                    new_transition.corner_IDs[i*(2*M)+2*j+1]=old_transition.corner_IDs[i*(2*M)+2*j+1];
                                }           
                            }
                        }

                        counter=0;
#pragma omp parallel
                        {
                            /* variables for managing the region */
                            std::array<abs_type,state_dim> lb;  /* lower-left corner */
                            std::array<abs_type,state_dim> ub;  /* upper-right corner */
                            std::array<abs_type,state_dim> no;  /* number of cells per dim */
                            std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */
                            /* second loop: fill pre array */
#pragma omp for schedule(dynamic, CHUNK_SIZE)
                            for(abs_type i=0; i<N; i++) {
                                /* loop over all inputs */
                                for(abs_type j=0; j<M; j++) {
                                    /* is x an element of the overflow symbols ? */
                                    if(new_transition.out_of_domain[i*M+j]) 
                                        continue;
                                    /* extract lower-left and upper-bound points */
                                    abs_type k_lb=new_transition.corner_IDs[i*(2*M)+2*j];
                                    abs_type k_ub=new_transition.corner_IDs[i*2*M+2*j+1];
                                    abs_type npost=1;

                                    /* cell idx to coordinates */
                                    for(int k=state_dim-1; k>=0; k--) {
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
                                        for(std::size_t l=0; l<state_dim; l++) 
                                            q+=(lb[l]+cc[l])*NN[l];
                                        cc[0]++;
                                        for(std::size_t l=0; l<state_dim-1; l++) {
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
                                progress(i,counter);
                            } // end for
                        } // end parallel

                        fill_pre(new_transition);

                        write_to_filebb(m_state_alphabet, m_input_alphabet, recomputed,"recomputation_lazy");

                    }//function closed

                bool region(abs_type q,state_type d_lb, state_type d_ub, state_type eta){
                    state_type x;
                    m_state_alphabet.itox(q,x);
                    bool belong=true;
                    for(std::size_t i=0; i<state_dim; i++){
                        if (d_lb[i]-eta[i]/2.0 > x[i] || x[i] > d_ub[i]+eta[i]/2.0)
                            belong = false;
                    }
                    return belong;
                }

               


                     template<class F2, class F3, class F4=decltype(params::avoid_abs)>
                    void simAbstraction(TransitionFunction& new_transition,
                            std::queue<abs_type>& diff,
                            const TransitionFunction& old_transition, 
                            F2& d_lb,
                            F2& d_ub,
                            F3& rs_post,
                            F4& avoid=params::avoid_abs){

                        
                        /* variables for managing the region */
                        std::array<abs_type,state_dim> lb;  /* lower-left corner */
                        std::array<abs_type,state_dim> ub;  /* upper-right corner */
                        std::array<abs_type,state_dim> no;  /* number of cells per dim */
                        std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */

                        /* init in transition_function the members no_pre, no_post, pre_ptr */ 
                        new_transition.init_infrastructure(N,M);
                        /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
                        /* make if a states recompute or not.*/
                        std::unique_ptr<std::atomic_flag[]> recomputed_mark(new std::atomic_flag[N*M]());
                        for (abs_type i = 0; i < N*M; i++) {
                            recomputed_mark[i].clear(); //ATOMIC_FLAG_INIT;
                        }
                        /*contain the states which need to recompute*/
                        std::vector<bool> recompute_queue(N,false); 
                        std::vector<bool> out_of_region(N, true); 
                        std::vector<bool> diff_done(N,false);

                        /*some variables for computing neighbours*/
                        std::vector<int> tmp;
                        std::vector<std::vector<int>>  cc_neighbours;

                        int num[]={-1,0,1};
                        /*init a queue which contains the states around (d_lb, d_ub), recompute transtions for queue, 
                         * if in ODE, x gets into the region(lb,ub), then enqueue the  neighbours. repeat this process untill 
                         * queue become empty. make all states in queue recomputed.
                         */
                        /*initialize recompute_queue with the bigger box contains around states and (d_lb, d_ub)*/
                        abs_type nNRegion=1;
                        TicToc time_;
                        time_.tic();
                        for (std::size_t i = 0; i < state_dim; ++i)
                        {
                            tmp.push_back(num[0]);
                        }
                        cc_neighbours.push_back(tmp);

                        for (int i = state_dim-1; i>=0; i--)
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
                        for(std::size_t k=0; k<state_dim; k++) {
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
                            for(std::size_t l=0; l<state_dim; l++)  {
                                q+=(lb[l]+cc[l])*m_state_alphabet.get_nn()[l];
                            }
                            cc[0]++;
                            for(std::size_t l=0; l<state_dim-1; l++) {
                                if(cc[l]==no[l]) {
                                    cc[l]=0;
                                    cc[l+1]++;
                                }
                            }
                            for(abs_type j = 0; j < M; j++)
                            {
                                bool already_done = recomputed_mark[q*M+j].test_and_set();
                                if (!already_done){
                                    recompute_queue[q]=true;
                                    
                                }
                            }
                            if(region(q,d_lb,d_ub,eta))
                                out_of_region[q]=false;  
                               
                            }     
                            time_.toc();
                        
                       std::cout<<"initial recomputing queue size:"<<recompute_queue.size()<<std::endl;


                        abs_type conn=0;
                        abs_type coun=0;
                        /*start big loop untill the recompute_queue become empty*/

                        /* for display purpose */
                        abs_type counter=0;

#pragma omp parallel shared(counter)
                        {
                            //std::cout<<"the number of threads:"<<omp_get_num_threads()<<std::endl;
                            state_type x;
                            input_type u;
                            ds_type y;
                            state_type r; 
                            /* variables for managing the region */
                            std::array<abs_type,state_dim> lb;  /* lower-left corner */
                            std::array<abs_type,state_dim> ub;  /* upper-right corner */
                            std::array<abs_type,state_dim> no;  /* number of cells per dim */
                            std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */

#pragma omp for schedule(dynamic, CHUNK_SIZE)
                            for(abs_type i=0; i<N; i++) {
                                if (!recompute_queue[i]) {
                                    continue;
                                }
                                
                               abs_type q=i;
                               conn++;
                                if(avoid(q)) {
                                    for(abs_type j=0; j<M; j++) {
                                        new_transition.out_of_domain[q*M+j]=true;
                                    }
                                    continue;
                                }
                                /* loop over all inputs */
                                for(abs_type j=0; j<M; j++) {
                                    /* is i an element of the avoid symbols ? */

                                    new_transition.out_of_domain[q*M+j]=false;
                                    /* get center x of cell */
                                    m_state_alphabet.itox(q,x);
                                    /* cell radius (including measurement errors) */
                                    for(std::size_t k=0; k<state_dim; k++){
                                        r[k]=eta[k]/2.0+m_z[k];
                                        y[k]=x[k];
                                        y[k+state_dim]=r[k];
                                    }

                                    /* current input */
                                    m_input_alphabet.itox(j,u);
                                    /* integrate system and radius growth bound */
                                    /* the result is stored in x and r */

                        
                                        rs_post(y,u); 
                                        abs_type npost=1;
                                        process_post(new_transition, q, j, x, y, r, npost, lb, ub, no, cc);

                                        if(new_transition.out_of_domain[q*M+j]){
                                            continue;
                                        }

                                        compute_post_indices(new_transition, q, j, npost, lb, no, cc);
                                        /* increment number of transitions by number of post */
                                        new_transition.m_no_post[q*M+j]=npost;
                                        
                                }
                                
                                /* print progress */
                                if(m_verbose) {
                                    if(counter==0)
                                        std::cout << "1st loop: ";
                                }
                                // std::cout << "Computed transition for state " << i << "\n" << std::endl;
                                progress(q,counter);
                            } // end for
                        } // end parallel
                        std::cout<<"size:"<<conn<<std::endl;
                        std::unique_ptr<bool[]> recomputed(new bool[N*M]());
                        for (abs_type i = 0; i < N*M; i++) {
                            recomputed[i] = recomputed_mark[i].test_and_set();
                        }
                        std::cout<<"total " <<coun<<std::endl;
                        /*copy from old transtions*/
                        for(abs_type i = 0; i < N; ++i )
                        { 
                            for (abs_type j = 0; j < M; ++j)
                            {
                                if(recomputed[i*M+j] && !diff_done[i] && !avoid(i)){
                                    diff.push(i);
                                    diff_done[i]=true;
                                }
                                //new_transition.out_of_domain[i*M+j]=standard_transition.out_of_domain[i*M+j];
                                else{
                                    new_transition.out_of_domain[i*M+j]=old_transition.out_of_domain[i*M+j];
                                    new_transition.m_no_post[i*M+j]=old_transition.m_no_post[i*M+j];
                                    new_transition.corner_IDs[i*(2*M)+2*j]=old_transition.corner_IDs[i*(2*M)+2*j];
                                    new_transition.corner_IDs[i*(2*M)+2*j+1]=old_transition.corner_IDs[i*(2*M)+2*j+1];
                                }           
                            }
                        }

                        counter=0;
#pragma omp parallel
                        {
                            /* variables for managing the region */
                            std::array<abs_type,state_dim> lb;  /* lower-left corner */
                            std::array<abs_type,state_dim> ub;  /* upper-right corner */
                            std::array<abs_type,state_dim> no;  /* number of cells per dim */
                            std::array<abs_type,state_dim> cc;  /* coordinate of current cell in the region */
                            /* second loop: fill pre array */
#pragma omp for schedule(dynamic, CHUNK_SIZE)
                            for(abs_type i=0; i<N; i++) {
                                /* loop over all inputs */
                                for(abs_type j=0; j<M; j++) {
                                    /* is x an element of the overflow symbols ? */
                                    if(new_transition.out_of_domain[i*M+j]) 
                                        continue;
                                    /* extract lower-left and upper-bound points */
                                    abs_type k_lb=new_transition.corner_IDs[i*(2*M)+2*j];
                                    abs_type k_ub=new_transition.corner_IDs[i*2*M+2*j+1];
                                    abs_type npost=1;

                                    /* cell idx to coordinates */
                                    for(int k=state_dim-1; k>=0; k--) {
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
                                        for(std::size_t l=0; l<state_dim; l++) 
                                            q+=(lb[l]+cc[l])*NN[l];
                                        cc[0]++;
                                        for(std::size_t l=0; l<state_dim-1; l++) {
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
                                progress(i,counter);
                            } // end for
                        } // end parallel

                        fill_pre(new_transition);

                        write_to_filebb(m_state_alphabet, m_input_alphabet, recomputed,"recomputation_lazy");

                    }//function closed

        };

} /* close namespace */
#endif /* ABSTRACTION_HH_ */
