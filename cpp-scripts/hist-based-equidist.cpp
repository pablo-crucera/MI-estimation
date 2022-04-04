#include <Rcpp.h>
#include <iostream>
#include <algorithm>    // std::sort
#include <math.h>

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
double histogram_based(NumericVector& x, NumericVector& y, int step = 3) {
    
    NumericVector X = clone(x);
    NumericVector Y = clone(y);

    if(X.size() != Y.size()) {
        std::cout << "Error: length of X and Y differ." << std::endl;
        return(0.0);
    }

    int N = X.size();

    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());

    NumericVector x_bins(N);
    NumericVector y_bins(N);

    NumericVector xlim(ceil(N % step) + 1);
    NumericVector ylim(ceil(N % step) + 1);

    double delta_x = x[1] - x[0];
    double delta_y = y[1] - y[0];

    xlim[0] = x[0];
    ylim[0] = y[0];

    for(int i = 1; i < N; i++) {
        x_bins[i] = i % step;
        y_bins[i] = i % step;
        if((i % step == 0) && (i != N - 1)) {
            xlim[i / step] = (x[i+1] + x[i]) / 2;
            ylim[i / step] = (y[i+1] + y[i]) / 2;
        }
    }
    
    IntegerVector idX(N);
    NumericVector idY(N);
    
    for(int i = 0; i < N; ++i) {
      idX[i] = i;
      idY[i] = i;
    }
    
    std::cout << x.size() << ',' << idX.size() << std::endl;

    NumericMatrix sortedX = sort(X, idX);
    NumericMatrix sortedY = sort(Y, idY);

    sortedX = sort(sortedX(_,1), x_bins);
    sortedY = sort(sortedY(_,1), y_bins);
    // 
    // x_bins = sortedX(_, 1);
    // y_bins = sortedY(_, 1);

    // xlim[N-1] = x[N-1];
    // ylim[N-1] = y[N-1];

    // double minX = min(X);
    // double maxX = max(X);
    // double minY = min(Y);
    // double maxY = max(Y);

    // k_x = ceil((maxX - minX) / delta_x);
    // k_y = ceil((maxY - minY) / delta_y);

    std::cout << x_bins << std::endl << y_bins << std::endl << std::endl;
    std::cout << X << std::endl << Y << std::endl << std::endl;

    // delta_x = (maxX - minX)/(double) k_x;
    // delta_x = (maxX - minX)/(double) k_x;


    // NumericVector hist_X(k_x);
    // NumericVector hist_Y(k_y);
    // NumericMatrix h(k_x, k_y);

    // NumericVector X_bins = (x - minX)/delta_x;
    // NumericVector Y_bins = (y - minY)/delta_y;

    // for(int i = 0; i < N; ++i) {
    //     int binX = (int) std::max(0.0, ceil(X_bins[i]) - 1);
    //     int binY = (int) std::max(0.0, ceil(Y_bins[i]) - 1);

    //     hist_X[binX] += 1;
    //     hist_Y[binY] += 1;
    //     h(binX, binY) += 1;
    // }

    // double MI = 0;

    // for(int i = 0; i < k_x; ++i) {
    //     for(int j = 0; j < k_y; ++j) {
    //         if (h(i,j) != 0) MI += (double) h(i,j)*log(h(i,j)*N/(hist_X[i]*hist_Y[j]));
    //     }
    // }

    // MI /= (double) N;

    return(0.0);
}
    
    // std::sort(X.begin(), X.end());
    // std::sort(Y.begin(), Y.end());

    // int bin = 0;
    // int cnt = 0;
    // double low_lim = X[0];
    // double upp_lim = low_lim + delta_x;

    // for(int i = 0; i < X.size(); ++i) {
    //     if(bin == k_x - 1) {
    //         hist_X[bin] = X.size() - i;
    //         i = X.size();
    //     }
    //     else if((X[i] >= low_lim) && (X[i] <= upp_lim)) cnt +=1;
    //     else {
    //         hist_X[bin] = cnt;
    //         bin += 1;
    //         cnt = 0;
    //         --i;
    //         low_lim = upp_lim;
    //         upp_lim += delta_x;
    //     }
    // }

    

    // bin = 0;
    // cnt = 0;
    // low_lim = Y[0];
    // upp_lim = low_lim + delta_y;

    // for(int j = 0; j < Y.size(); ++j) {
    //     if(bin == k_x - 1) {
    //         hist_Y[bin] = Y.size() - j;
    //         j = Y.size();
    //     }
    //     else if((Y[j] >= low_lim) && (Y[j] <= upp_lim)) cnt +=1;
    //     else {
    //         hist_Y[bin] = cnt;
    //         bin += 1;
    //         cnt = 0;
    //         --j;
    //         low_lim = upp_lim;
    //         upp_lim += delta_y;
    //     }
    // }

    // std::cout << "hist_X: " << hist_X << std::endl;
    // std::cout << "hist_Y: " << hist_Y << std::endl;