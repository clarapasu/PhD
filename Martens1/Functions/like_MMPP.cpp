// MMPP likelihood for continuous time mark recapture
// based on code by Paul Blackwell
// C++ by David L Miller, GPLv2+
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() () {
  
  PARAMETER_VECTOR(theta);
  
  // actual data
  DATA_IVECTOR(Us);
  DATA_VECTOR(t);
  DATA_IVECTOR(s);
  DATA_VECTOR(f);
  
  DATA_SCALAR(t0);
  DATA_SCALAR(Time);
  DATA_INTEGER(n_states);
  DATA_INTEGER(n_indiv);
  DATA_INTEGER(n_obs);
  
  DATA_IVECTOR(Q_template);
  DATA_VECTOR(Q_fixed);
  DATA_IVECTOR(lambda_template);
  DATA_VECTOR(lambda_fixed);
  
  Type ll = 0.0;
  vector<Type> etheta = exp(theta);
  
  
  vector<Type> lambda(n_states);
  matrix<Type> Lambda(n_states, n_states);
  matrix<Type> Q(n_states, n_states);
  matrix<Type> A(n_states, n_states);
  
  // build Q
  int Qit = 0;
  for(int j=0; j<n_states; j++){
    for(int k=0; k<n_states; k++){
      if(Q_fixed(Qit) == 1){
        Q(j,k) = 0.0;
      }else{
        Q(j,k) = etheta(Q_template(Qit)-1);
      }
      Qit++;
    }
  }
  // set the diagonal elements
  Q.diagonal() << -1.0 * Q.rowwise().sum();
  
  // build lambda, this is less horrible
  for(int j=0; j<n_states; j++){
    if(lambda_fixed(j) == 1){
      lambda(j) = 0.0;
    }else{
      lambda(j) = etheta(lambda_template(j)-1);
    }
  }
  
  // create (diagonal) Lambda matrix
  Lambda.setZero();
  Lambda.diagonal() << lambda;
  
  matrix<Type>G = Q-Lambda;
  
  int istart = 0;
  
  Type ill;
  // loop over individuals
  for(int i=0; i<n_indiv; i++){
    
    // timeseries length for this individual
    int U = Us(i);
    
    // get data for this individual
    vector<Type> ti = t.segment(istart, U);
    vector<int> si = s.segment(istart, U);
    istart = istart + U;
    
    if (U==0){
      A = atomic::expm(matrix<Type>((Time-t0)*G));
      ill = (f.matrix().transpose()*A).sum();
    }else{
      vector<Type> ts(ti.size()+2);
      ts << t0, ti, Time;
      
      vector<Type> tau = diff(ts);
      vector<Type> w(U+1);
      
      // First interval - starts with possibly unknown state
      A = atomic::expm(matrix<Type>(tau(0)*G));
      
      vector<Type> omega1 = A.col(si(0));
      w(0) = (omega1.array()*f.array()).sum();
      w(0) = w(0) * lambda(si(0));
      
      if(U>1){
        for(int u=1; u<U; u++){
          A = atomic::expm(matrix<Type>(tau(u)*G));
          w(u) = A(si(u-1), si(u)) * lambda(si(u));
        }
      }
      
      // Final interval - end state unknown
      A = atomic::expm(matrix<Type>(tau(U)*G));
      vector<Type> omega_Uplus1 = A.row(si(U-1));
      w(U) = sum(omega_Uplus1);
      ill = w.prod();
    }
    
    // add to likelihood
    ll -= log(ill);
  } 
  
  
  vector<Type> ones(n_states);
  ones.setOnes();
  A = atomic::expm(matrix<Type>(Time * G));
  Type p_unobs = (f.matrix().transpose() * A * ones.matrix())(0,0);
  Type p_obs = 1.0 - p_unobs;
  
  ll += n_obs * log(p_obs);
  
  // want the SEs of exp(theta), this avoids us having to delta method
  ADREPORT(etheta);
  
  return(ll);
} 
