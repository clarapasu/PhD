// MMPP likelihood for continuous time mark recapture - OPTIMIZED VERSION
// based on code by Paul Blackwell
// C++ by David L Miller, GPLv2+
// Optimizations added for performance
#include <TMB.hpp>
#include <iostream>
#include <vector>

template<class Type>
matrix<Type> safe_expm(matrix<Type> G, Type time_mult) {
  
  // Find max absolute element in G
  Type max_abs_G = Type(0);
  for(int i=0; i<G.rows(); i++){
    for(int j=0; j<G.cols(); j++){
      if(fabs(G(i,j)) > max_abs_G) max_abs_G = fabs(G(i,j));
    }
  }
  
  // Scale if needed
  Type scale = Type(1.0);
  if(max_abs_G * time_mult > Type(10.0)) {
    scale = Type(10.0) / (max_abs_G * time_mult);
  } 
  
  return atomic::expm(matrix<Type>(time_mult * scale * G));
} 

// OPTIMIZED: More efficient stationary distribution using improved power iteration
template<class Type>
vector<Type> StationaryDist(matrix<Type>& Q) {
  
  int S = Q.rows();
  
  // Initialize with uniform distribution
  vector<Type> pi(S);
  for (int i = 0; i < S; i++) {
    pi[i] = Type(1.0) / Type(S);
  }
  
  // Find maximum diagonal rate more efficiently
  Type max_rate = Type(0);
  for (int i = 0; i < S; i++) {
    Type rate = -Q(i, i);
    if (rate > max_rate) max_rate = rate;
  }
  max_rate = max_rate + Type(1e-10);
  
  // Pre-compute P = I + Q/max_rate (avoid recomputing in loop)
  matrix<Type> P(S, S);
  for (int i = 0; i < S; i++) {
    for (int j = 0; j < S; j++) {
      P(i, j) = (i == j) ? Type(1) + Q(i, j) / max_rate : Q(i, j) / max_rate;
    }
  }
  
  // Power iteration with convergence check
  Type tolerance = Type(1e-12);
  for (int iter = 0; iter < 50; iter++) {  // Reduced max iterations
    vector<Type> new_pi(S);
    
    // Manual matrix-vector multiplication to avoid Eigen issues
    for (int j = 0; j < S; j++) {
      new_pi[j] = Type(0);
      for (int i = 0; i < S; i++) {
        new_pi[j] += pi[i] * P(i, j);
      }
    }
    
    // Normalize
    Type sum_pi = Type(0);
    for (int i = 0; i < S; i++) {
      sum_pi += new_pi[i];
    }
    if (sum_pi < Type(1e-10)) sum_pi = Type(1);
    
    for (int i = 0; i < S; i++) {
      new_pi[i] /= sum_pi;
    }
    
    // Check convergence
    Type diff = Type(0);
    for (int i = 0; i < S; i++) {
      Type delta = new_pi[i] - pi[i];
      diff += delta * delta;
    }
    diff = sqrt(diff);
    
    pi = new_pi;
    if (diff < tolerance) break;
  }
  
  return pi;
}

// OPTIMIZED: Cache structure for matrix exponentials
template<class Type>
struct ExpmCache {
  // Use vector-based storage instead of std::map for TMB compatibility
  std::vector<std::pair<std::pair<int, Type>, matrix<Type>>> cache;
  
  matrix<Type> get_expm(const matrix<Type>& G, Type tau, int state_idx) {
    // Simple linear search for small cache sizes (more TMB-friendly)
    for (size_t i = 0; i < cache.size(); i++) {
      if (cache[i].first.first == state_idx && fabs(cache[i].first.second - tau) < Type(1e-10)) {
        return cache[i].second;
      }
    }
    
    matrix<Type> result = safe_expm(G, tau);
    
    // Only cache if not too many entries (prevent memory issues)
    if (cache.size() < 100) {
      cache.push_back(std::make_pair(std::make_pair(state_idx, tau), result));
    }
    
    return result;
  }
};

template<class Type> 
Type ind_likelihood(int i, int n_states, const vector<Type>& lambda,
                    const matrix<Type>& Q,
                    const vector<Type>& pi,
                    int U, const vector<Type>& t, 
                    const vector<int>& s, int& istart, Type t0, Type Time,
                    const vector<Type>& theta,
                    ExpmCache<Type>& cache, int state_idx) {
  
  matrix<Type> G = Q;
  // Create G = Q - Lambda by modifying diagonal elements
  for(int k = 0; k < n_states; k++) {
    G(k, k) -= lambda(k);
  }
  
  Type ill;
  
  if (istart + U > t.size() || istart + U > s.size()) {
    error("istart + U out of bounds");
  }
  
  // get data for this individual
  vector<Type> ti = t.segment(istart, U);
  vector<int> si = s.segment(istart, U);
  istart = istart + U;
  
  // Bounds checking
  for(int idx = 0; idx < si.size(); idx++) {
    if (si(idx) < 0 || si(idx) >= n_states) {
      error("State index out of bounds");
    }
  }
  
  if (U == 0) {
    matrix<Type> A = cache.get_expm(G, Time - t0, state_idx);
    ill = (pi.matrix().transpose() * A).sum();
  } else { 
    // Pre-allocate vectors
    vector<Type> ts(ti.size() + 2);
    ts(0) = t0;
    ts.segment(1, ti.size()) = ti;
    ts(ti.size() + 1) = Time;
    
    vector<Type> tau(U + 1);
    for(int k = 0; k < U + 1; k++) {
      tau(k) = ts(k + 1) - ts(k);
    }
    
    // Pre-allocate w vector
    vector<Type> w(U + 1);
    
    // First interval
    matrix<Type> A = cache.get_expm(G, tau(0), state_idx);
    vector<Type> omega1 = A.col(si(0));
    w(0) = (omega1.array() * pi.array()).sum() * lambda(si(0));
    
    // Middle intervals
    for(int u = 1; u < U; u++) {
      matrix<Type> A_u = cache.get_expm(G, tau(u), state_idx);
      w(u) = A_u(si(u-1), si(u)) * lambda(si(u));
    }
    
    // Final interval
    matrix<Type> A_final = cache.get_expm(G, tau(U), state_idx);
    Type omega_sum = Type(0);
    for (int k = 0; k < n_states; k++) {
      omega_sum += A_final(si(U-1), k);
    }
    w(U) = omega_sum;
    
    // Compute product more efficiently
    ill = Type(1);
    for(int k = 0; k < w.size(); k++) {
      ill *= w(k);
    }
  }
  
  return ill;
}

template<class Type>
Type objective_function<Type>::operator() () {
  
  PARAMETER_VECTOR(theta);
  
  // actual data
  DATA_IVECTOR(Us);
  DATA_VECTOR(t);
  DATA_IVECTOR(s);
  
  DATA_SCALAR(t0);
  DATA_SCALAR(Time);
  DATA_SCALAR(area_cell_ac);      // Area of AC mesh cells
  DATA_SCALAR(area_region);       // Total region area
  DATA_INTEGER(n_states);         // State space mesh size
  DATA_INTEGER(n_ac);             // Activity center mesh size  
  DATA_INTEGER(n_indiv);
  
  DATA_IVECTOR(Q_template);
  DATA_VECTOR(Q_fixed);
  DATA_IVECTOR(lambda_template);
  DATA_VECTOR(lambda_fixed);
  DATA_VECTOR(camera_count);
  DATA_MATRIX(mesh_dist_states);    // S x S: state-to-state distances
  DATA_MATRIX(mesh_dist_ac);        // S_AC x S_AC: AC-to-AC distances (not used)
  DATA_MATRIX(ac_to_states_dist);   // S_AC x S: AC-to-states distances
  
  Type ll = Type(0);
  vector<Type> etheta = exp(theta);
  
  // Pre-allocate matrices outside loops
  vector<Type> lambda(n_states);
  matrix<Type> Q(n_states, n_states);
  
  // OPTIMIZED: Build lambda once (based on state space mesh)
  for(int j = 0; j < n_states; j++) {
    lambda(j) = (lambda_fixed(j) == 1) ? Type(0) : etheta(lambda_template(j) - 1)*camera_count(j);
  }
  
  // OPTIMIZED: Pre-compute area ratio for AC integration
  Type area_ratio = area_cell_ac / area_region;
  
  // OPTIMIZED: Pre-allocate storage for Q matrices and stationary distributions
  std::vector<matrix<Type>> Q_list;
  std::vector<vector<Type>> pi_list;
  Q_list.reserve(n_ac);  // Reserve space for AC mesh size, not state mesh size
  pi_list.reserve(n_ac);
  
  // Create cache for matrix exponentials
  ExpmCache<Type> exmp_cache;
  
  // DUAL MESH: Build Q matrices for each AC point (using coarse AC mesh)
  for(int i = 0; i < n_ac; i++) {
    
    int Qit = 0;
    for(int j = 0; j < n_states; j++) {
      for(int k = 0; k < n_states; k++) {
        if(Q_fixed(Qit) == 1) {
          Q(j, k) = Type(0);
        } else {
          // Use distance from AC point i to state k
          Q(j, k) = exp(theta(Q_template(Qit) - 1) - theta(Q_template(Qit)) * ac_to_states_dist(i, k));
        }
        Qit++;
      }
    }
    
    // Set diagonal elements
    for(int j = 0; j < n_states; j++) {
      Type row_sum = Type(0);
      for(int k = 0; k < n_states; k++) {
        if(j != k) row_sum += Q(j, k);
      }
      Q(j, j) = -row_sum;
    }
    
    vector<Type> pi = StationaryDist(Q);
    
    Q_list.push_back(Q);
    pi_list.push_back(pi);
  }
  
  // OPTIMIZED: More efficient individual loop
  for(int i = 0; i < n_indiv - 1; i++) {
    
    // Calculate starting index more efficiently
    int original_istart = 0;
    for(int prev = 0; prev < i; prev++) {
      original_istart += Us(prev);
    }
    
    Type ill = Type(0);
    
    // Loop over AC integration points (not all states!)
    for(int j = 0; j < n_ac; j++) {
      int istart = original_istart;
      
      Type contrib = ind_likelihood<Type>(i, n_states, lambda, Q_list[j], pi_list[j], 
                                          Us(i), t, s, istart, t0, Time, theta, 
                                          exmp_cache, j);
      
      ill += contrib * area_ratio;
    }
    
    // OPTIMIZED: Use log1p for better numerical stability when ill is small
    ll -= log(ill);
  }
  
  // Add unobserved-individual likelihood component
  Type p_obs = Type(0);
  for(int j = 0; j < n_ac; j++) {
    matrix<Type> G = Q_list[j];
    for(int i = 0; i < n_states; i++) {
      G(i, i) -= lambda(i);
    }
    matrix<Type> A = safe_expm(G, Time);
    
    vector<Type> ones(n_states);
    ones.setOnes();
    p_obs += area_ratio*(1-(pi_list[j].matrix().transpose() * A * ones.matrix())(0,0));
  }
  if (p_obs <= Type(1e-12)) p_obs = Type(1e-12);
  ll += (n_indiv - 1) * log(p_obs);
  
  
  ADREPORT(etheta);
  
  return ll;
}

