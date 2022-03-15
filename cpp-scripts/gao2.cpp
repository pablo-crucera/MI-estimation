#include <Rcpp.h>
#include <iostream>
#include <algorithm>    // std::sort
#include <math.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/digamma.hpp>
// #include <dplyr.h>
// // [[Rcpp::depends(dplyr, plogr, bindrcpp)]]

using namespace Rcpp ;


// [[Rcpp::export]]
NumericMatrix merge(NumericMatrix M1, NumericMatrix M2){
  NumericVector A = M1(_, 0);
  NumericVector B = M2(_, 0);

  NumericVector idx_A = M1(_, 1);
  NumericVector idx_B = M2(_, 1);
  
  int N = A.size();
  int M = B.size();
  int i = 0;
  int j = 0;
  
  NumericVector out = NumericVector(N + M);
  NumericVector idx_out = NumericVector(N + M);
  
  int cnt = 0;
  
  while((i != N) || (j != M)) {
    if(i == N) {
      for(int k = j; k < M; k++) {
        out[k + N] = B[k];
        idx_out[k + N] = idx_B[k];
      }
      j = M;
    } 
    else if(j == M) {
      for(int k = i; k < N; k++) {
        out[k + M] = A[k];
        idx_out[k + M] = idx_A[k];
      }
      i = N;
    } 
    else if(A[i] >= B[j]){ 
      out[cnt] = B[j];
      idx_out[cnt] = idx_B[j];
      j++;
    } 
    else {
      out[cnt] = A[i];
      idx_out[cnt] = idx_A[i];
      i++;
    }
    ++cnt;
  }
  
  NumericMatrix outM( N + M , 2);
  NumericMatrix::Column col1 = outM(_, 0);
  NumericMatrix::Column col2 = outM(_, 1);

  col1 = out;
  col2 = idx_out;
  return(outM);
}

// [[Rcpp::export]]
NumericMatrix sort(NumericVector v, NumericVector idx){
  NumericMatrix outM(v.size(), 2);
  NumericMatrix::Column col1 = outM(_, 0);
  NumericMatrix::Column col2 = outM(_, 1);
  
  if(v.size() == 1){
    col1 = v;
    col2 = idx;
    return(outM);
  } 
  else if(v.size() == 2) {
    if(v[1] >= v[0]) {
      col1 = v;
      col2 = idx;
      return(outM);
    } 
    else {
      col1 = NumericVector::create(v[1], v[0]);
      col2 = NumericVector::create(idx[1], idx[0]);
      return(outM);
    } 
  } else if(v.size() == 0) return(as<NumericMatrix>(0));
  else return(merge(
    sort(
      head(v, ceil((float) v.size() / (float) 2 - 0.01)),
      head(idx, ceil((float) v.size() / (float) 2 - 0.01))
    ),
    sort(
      tail(v, floor((float) v.size() / 2.0 + 0.01)),
      tail(idx, floor((float) v.size() / 2.0 + 0.01))
    )
  ));
}

// [[Rcpp::export]]

NumericVector mapper(NumericVector v, NumericVector idx) {
  NumericVector out(v.size());

  for(int i = 0; i < v.size(); i++) out[i] = v[idx[i]];

  return(out);
}



// [[Rcpp::export]]

double gao2(NumericVector X, NumericVector Y, int k = 5, int max_neighbors = 30) {
  if(X.size() != Y.size()) {
    std::cout << "Error: sizes of X and Y must be the same"  << std::endl;
    return(0);
  } else {
    int N = X.size();
    NumericMatrix dist_X(N, N);
    NumericMatrix dist_Y(N, N);
    NumericMatrix dist(N, N);
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        dist_X(i,j) = std::abs(X(i) - X(j));
        dist_Y(i,j) = std::abs(Y(i) - Y(j));
        dist(i,j) = std::max(dist_X(i,j), dist_Y(i,j));
      }
    }
    
    // std::cout << "-----------------------------------------" << std::endl << "dist" << std::endl << dist
              // << std::endl << "-----------------------------------------" << std::endl;

    double xi = 0;
    
    for(int i = 0; i < N; i++) {
      NumericMatrix sorted_dist_i(N, 3);
      NumericMatrix::Column col1 = sorted_dist_i(_, 0);
      NumericMatrix::Column col2 = sorted_dist_i(_, 1);
      NumericMatrix::Column col3 = sorted_dist_i(_, 2);
      
      NumericVector forprinting = dist(_,i);
      

      NumericVector idx(0);
      for(int i = 0; i < N; i++) idx.insert(i,i);
      NumericMatrix sorted_dist_idx = sort(dist(_,i), idx);
      NumericVector sorted_idx = sorted_dist_idx(_, 1);
      // std::cout<< "-----------------------------------------" << std::endl << "sorted_dist_idx" << std::endl << sorted_dist_idx
      //         << std::endl << "-----------------------------------------" << std::endl;
      col1 = sorted_dist_idx(_, 0);
      col2 = mapper(dist_X(_, i), sorted_idx);
      col3 = mapper(dist_Y(_, i), sorted_idx);
      
      // std::cout<< "-----------------------------------------" << std::endl << "sorted_dist_i" << std::endl << sorted_dist_i
      //         << std::endl << "-----------------------------------------" << std::endl;
      
      double rho_i = sorted_dist_i(k, 0);
      double rho_xi = sorted_dist_i(k, 1);
      double rho_yi = sorted_dist_i(k, 2);
      int k_i = 0;
      
      if (rho_i != 0) k_i = k;
      else {
        for(int j = 1; sorted_dist_i[j] == 0; j++) k_i++;
      }
      
      NumericVector sorted_dist_xi = dist_X(_, i);
      std::sort(sorted_dist_xi.begin(), sorted_dist_xi.begin()+N);
      NumericVector sorted_dist_yi = dist_Y(_, i);
      std::sort(sorted_dist_yi.begin(), sorted_dist_yi.begin() + N);

      int n_xi = 0; int n_yi = 0;
      
      for(int j = 1; sorted_dist_xi[j] <= rho_xi && n_xi <= max_neighbors; j++) n_xi++;
      for(int j = 1; sorted_dist_yi[j] <= rho_yi && n_yi <= max_neighbors; j++) n_yi++;
      //std::cout<< "(N, n_xi, n_yi, k_i): (" << N << "," << n_xi << "," << n_yi << "," << k_i << ")" << std::endl;
      xi += boost::math::digamma(k_i) + log(N) - log(n_xi + 1) - log(n_yi + 1);
      //std::cout<< "Xi = (" << boost::math::digamma(k_i) << ", " << log(N) << ", " << log(n_xi + 1) << ", " << log(n_yi + 1) << ")" << std::endl << std::endl;
    }
    
    return(xi / (double) N);
  }
}