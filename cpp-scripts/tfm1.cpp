#include <Rcpp.h>
#include <iostream>
#include <algorithm>    // std::sort
#include <math.h>

// [[Rcpp::depends(BH)]]
#include<boost/math/special_functions/digamma.hpp>
using namespace Rcpp;

// [[Rcpp::export]]


double gao1(NumericVector X, NumericVector Y, int k = 5, int max_neighbors = 30) {
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

        double xi = 0;

        for(int i = 0; i < N; i++) {
            NumericVector sorted_dist_i = dist(_, i);
            std::sort(sorted_dist_i.begin(), sorted_dist_i.begin()+N);
            NumericVector sorted_dist_xi = dist_X(_, i);
            std::sort(sorted_dist_xi.begin(), sorted_dist_xi.begin()+N);
            NumericVector sorted_dist_yi = dist_Y(_, i);
            std::sort(sorted_dist_yi.begin(), sorted_dist_yi.begin() + N);
            
            //std::cout<< "-----------------------------------------" << std::endl << sorted_dist_i
            //         << std::endl << sorted_dist_xi << std::endl << sorted_dist_yi << std::endl << "-----------------------" << std::endl;
            
            double rho_i = sorted_dist_i[k];
            int k_i = 0;

            if (rho_i != 0) k_i = k;
            else {
                for(int j = 1; sorted_dist_i[j] == 0; j++) k_i++;
            }

            int n_xi = 0; int n_yi = 0;

            for(int j = 1; sorted_dist_xi[j] <= rho_i && n_xi <= max_neighbors; j++) n_xi++;
            for(int j = 1; sorted_dist_yi[j] <= rho_i && n_yi <= max_neighbors; j++) n_yi++;
            //std::cout<< "(N, n_xi, n_yi, k_i): (" << N << "," << n_xi << "," << n_yi << "," << k_i << ")" << std::endl;
            xi += boost::math::digamma(k_i) + log(N) - log(n_xi + 1) - log(n_yi + 1);
            //std::cout<< "Xi = (" << boost::math::digamma(k_i) << ", " << log(N) << ", " << log(n_xi + 1) << ", " << log(n_yi + 1) << ")" << std::endl << std::endl;
        }

        return(xi / (double) N);
        
    }
}

