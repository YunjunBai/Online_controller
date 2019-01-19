
#ifndef SOLVER_HH_
#define SOLVER_HH_

#include <iostream>
#include <climits>
#include <stdexcept>
#include <queue>
#include <memory>
#include <utility>

#include "UniformGrid.hh"
#include "TransitionFunction.hh"

#include "WinningDomain.hh"

namespace scots{


/******************solver for safety game***************************/
template<class F>
WinningDomain static_safety_game(const TransitionFunction& trans_function, F& safe) {
  /* size of state alphabet */
  abs_type N=trans_function.m_no_states;
  /* size of input alphabet */
  abs_type M=trans_function.m_no_inputs;
  /* used to encode that a state is not in the winning domain */
  abs_type loosing = std::numeric_limits<abs_type>::max();

  /* valid_inputs
   * boolean array of size N*M
   * valid_inputs[i*M+j]=true iff input j
   * is a valid input at state i */
  std::vector<bool> valid_inputs(N*M,false); 
  /* no_input: keep track of the number of valid inputs.
   * If no_val_in[i]==0, then state i is not winning and 
   * no_val_in[i] is set to loosing (see also WinningDomain) */
  std::vector<abs_type> no_val_in(N,0);

  std::vector<double> value;
   value.resize(N,std::numeric_limits<double>::infinity());
   std::vector<double> value_temp(N,0);
  std::unique_ptr<double[]>  edge_val(new double[N*M]);
  /* initialization */
  std::queue<abs_type> fifo;
  for(abs_type i=0; i<N; i++) {
    if(safe(i)) {
      for(abs_type j=0; j<M; j++) {
        if(trans_function.m_no_post[i*M+j]) {
          valid_inputs[i*M+j]=true;
          no_val_in[i]++;
        }
      }
    }
    if(!no_val_in[i]) {
      fifo.push(i);
      value[i]=0;
      /* mark no_val_in[i]=loosing to indicate that state i ist not winning */
      no_val_in[i]=loosing;
    }
  }

  while(!fifo.empty()) {
    abs_type k=fifo.front();
    fifo.pop();
    /* loop over all inputs */
    for(abs_type j=0; j<M; j++) {
      /* loop over all pre states of (k,j) */
      for(abs_ptr_type p=0; p<trans_function.m_no_pre[k*M+j]; p++) {
        /* (i,j,k) is a transition */
        abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[k*M+j]+p];

        edge_val[i*M+j]=(edge_val[i*M+j]>=1+value[k] ? 1+value[k] : edge_val[i*M+j]);
        
        value_temp[i]=(edge_val[i*M+j]+1>=value_temp[i] ? edge_val[i*M+j]+1 : value_temp[i]);
        /* check if input j at state i is considered safe */
        if(valid_inputs[i*M+j]) {
          /* set source states i with label j as unsafe pair */
          valid_inputs[i*M+j]=false;
          no_val_in[i]--;
        }
        /* add to unsafe set if state i has no more valid inputs */
        if(!no_val_in[i] && value[i]!=value_temp[i]) {
          fifo.push(i);
          value[i]=value_temp[i];
          /* mark no_val_in[i]=loosing to indicate that state i ist not winning */
          no_val_in[i]=loosing;
        }
      }

    }
  }
  
  return WinningDomain(N,M,std::move(no_val_in),std::move(value),std::move(valid_inputs),loosing);
}


WinningDomain online_safety_game(const TransitionFunction& trans_function, 
                                 std::queue<abs_type> online_queue, 
                                 WinningDomain wd) {
  /* size of state alphabet */
  abs_type N=trans_function.m_no_states;
  /* size of input alphabet */
  abs_type M=trans_function.m_no_inputs;
  /* used to encode that a state is not in the winning domain */
  abs_type loosing = std::numeric_limits<abs_type>::max();

   /* initialization */
  std::vector<bool> valid_inputs=wd.get_m_inputs(); 
  
  std::vector<abs_type> no_val_in=wd.get_m_winning_domain();

  std::vector<double> value=wd.get_value();

  std::vector<double> value_temp(N,0);
  std::unique_ptr<double[]>  edge_val(new double[N*M]);
 
  std::queue<abs_type> fifo=online_queue;
  

  while(!fifo.empty()) {
    abs_type k=fifo.front();
    fifo.pop();
    /* loop over all inputs */
    for(abs_type j=0; j<M; j++) {
      /* loop over all pre states of (k,j) */
      for(abs_ptr_type p=0; p<trans_function.m_no_pre[k*M+j]; p++) {
        /* (i,j,k) is a transition */
        abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[k*M+j]+p];

        edge_val[i*M+j]=(edge_val[i*M+j]>=1+value[k] ? 1+value[k] : edge_val[i*M+j]);
        
        value_temp[i]=(edge_val[i*M+j]+1>=value_temp[i] ? edge_val[i*M+j]+1: value_temp[i]);
        /* check if input j at state i is considered safe */
        if(valid_inputs[i*M+j]) {
          /* set source states i with label j as unsafe pair */
          valid_inputs[i*M+j]=false;
          no_val_in[i]--;
        }
        /* add to unsafe set if state i has no more valid inputs */
        if(!no_val_in[i] && value[i]!=value_temp[i]) {
          fifo.push(i);
          value[i]=value_temp[i];
          /* mark no_val_in[i]=loosing to indicate that state i ist not winning */
          no_val_in[i]=loosing;
        }
        
      }
    }
  }
  return WinningDomain(N,M,std::move(no_val_in),std::move(value),std::move(valid_inputs),loosing);
}




/********************solver for reachability game*********************/
template<class F1, class F2=decltype(params::avoid)>
WinningDomain static_reachability_game(const TransitionFunction& trans_function,
                                      F1& target, 
                                      F2& avoid = params::avoid) {
  /* size of state alphabet */
  abs_type N=trans_function.m_no_states;
  /* size of input alphabet */
  abs_type M=trans_function.m_no_inputs;

  /* used to encode that a state is not in the winning domain */
  abs_type loosing = std::numeric_limits<abs_type>::max();
  if(M > loosing-1) {
    throw std::runtime_error("scots::solve_reachability_game: Number of inputs exceeds maximum supported value");
  }
  /* win_domain[i] = j
   * contains the input j associated with state i
   *
   * j = loosing if the target is not reachable from i 
   *
   * initialize all states to loosing (win_domain[i]=loosing) */
  std::vector<abs_type> win_domain(N,loosing); 
  /* initialize value */
  std::vector<double> value;
  value.resize(N,std::numeric_limits<double>::infinity());
  /* keep track of the number of processed post */
  std::unique_ptr<abs_type[]> K(new abs_type[N*M]);
  /* keep track of the values (corresponds to M in Alg.2)*/
  std::unique_ptr<double[]>  edge_val(new double[N*M]);

  /* init fifo */
  std::queue<abs_type> fifo;
  for(abs_type i=0; i<N; i++) {
    if(target(i) && !avoid(i)) {
      win_domain[i]=loosing;
      /* value is zero */
      value[i]=0;
      /* states in the target are added to the fifo */
      fifo.push(i);
    }
    for(abs_type j=0; j<M; j++) {
      edge_val[i*M+j]=0;
      K[i*M+j]=trans_function.m_no_post[i*M+j];
    }
  }

  /* main loop */
  while(!fifo.empty()) {
    /* get state to be processed */
    abs_type q=fifo.front();
    fifo.pop();
    /* loop over each input */
    for(abs_type j=0; j<M; j++) {
      /* loop over pre's associated with this input */
      for(abs_ptr_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
        abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[q*M+j]+v];
        if(avoid(i))
          continue;
        /* (i,j,q) is a transition */
        /* update the number of processed posts */
        K[i*M+j]--;

        /* update the max value of processed posts */
        edge_val[i*M+j]=(edge_val[i*M+j]>=1+value[q] ? edge_val[i*M+j] : 1+value[q]);
      
        /* check if for node i and input j all posts are processed */
        if(!K[i*M+j] && value[i]>edge_val[i*M+j]) {

          fifo.push(i);
          value[i]=edge_val[i*M+j]; 
          win_domain[i]=j;

        }
      }  /* end loop over all pres of state i under input j */
    }  /* end loop over all input j */
  }  /* fifo is empty */
  
  return WinningDomain(N,M,std::move(win_domain),std::move(value),std::vector<bool>{},loosing);
}



template<class F1=decltype(params::avoid)>
WinningDomain online_reachability_game(const TransitionFunction& trans_function,
                                      std::queue<abs_type> online_queue, 
                                      F1& avoid,
                                      WinningDomain wd
                                      ) {
  /* size of state alphabet */
  abs_type N=trans_function.m_no_states;
  /* size of input alphabet */
  abs_type M=trans_function.m_no_inputs;

  /* used to encode that a state is not in the winning domain */
  abs_type loosing = std::numeric_limits<abs_type>::max();
  if(M > loosing-1) {
    throw std::runtime_error("scots::solve_reachability_game: Number of inputs exceeds maximum supported value");
  }
  /* win_domain[i] = j
   * contains the input j associated with state i
   *
   * j = loosing if the target is not reachable from i 
   *
   * initialize all states to old win_domain (win_domain[i]) */
  std::vector<abs_type> win_domain=wd.get_m_winning_domain(); // todo 
  /* initialize value with old value */
  std::vector<double> value=wd.get_value(); //todo
  
  /* keep track of the number of processed post */
  std::unique_ptr<abs_type[]> K(new abs_type[N*M]); //todo
  /* keep track of the values (corresponds to M in Alg.2)*/
  std::unique_ptr<double[]>  edge_val(new double[N*M]); //todo  i am not sure

  /* init fifo */
  std::queue<abs_type> fifo=online_queue;

  for(abs_type i=0; i<N; i++) {  
    for(abs_type j=0; j<M; j++) {
      edge_val[i*M+j]=0; //todo
      K[i*M+j]=trans_function.m_no_post[i*M+j];
    }
  }

  /* main loop */
  while(!fifo.empty()) {
    /* get state to be processed */
    abs_type q=fifo.front();
    fifo.pop();
    /* loop over each input */
    for(abs_type j=0; j<M; j++) {
      /* loop over pre's associated with this input */
      for(abs_ptr_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
        abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[q*M+j]+v];
        if(avoid(i))
          continue;
        /* (i,j,q) is a transition */
        /* update the number of processed posts */
        K[i*M+j]--;
        /* update the max value of processed posts */
        edge_val[i*M+j]=(edge_val[i*M+j]>=1+value[q] ? edge_val[i*M+j] : 1+value[q]);
        /* check if for node i and input j all posts are processed */
        if(!K[i*M+j] && value[i]>edge_val[i*M+j]) {

          fifo.push(i);
          value[i]=edge_val[i*M+j]; 
          win_domain[i]=j;
        }
      }  /* end loop over all pres of state i under input j */
    }  /* end loop over all input j */
  }  /* fifo is empty */

  /* if the default value function was used, free the memory of the static object*/
  /*if(value == scots::params::value){
      value.clear();
      value.shrink_to_fit();
  }*/

  return WinningDomain(N,M,std::move(win_domain),std::move(value),std::vector<bool>{},loosing);
}

} /*namespace closed*/

#endif