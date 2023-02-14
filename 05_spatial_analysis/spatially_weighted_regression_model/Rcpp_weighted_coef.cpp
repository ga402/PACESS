#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * First function to generate the coef
 *
 */
// [[Rcpp::export]]
arma::mat rcpp_coefficient(const arma::mat X, const arma::vec w, arma::colvec Y) {
  arma::colvec diagonal(w.size());
  for (int i = 0; i < w.size(); i++){
    diagonal(i) = w[i];
  }
  return (inv(X.t()*diagmat(diagonal)*X)*X.t()*diagmat(diagonal)*Y).t();
}

/*
 * Run a loop
 */

// [[Rcpp::export]]
Rcpp::List rcpp_calculate_coefficient(const arma::mat kernel, const arma::mat X, arma::colvec Y){
  int ncols = X.n_cols;
  int nrows = kernel.n_rows;
  // the output matrix for the results  
  arma::mat coef_result(nrows, ncols);
  arma::mat result(1, nrows);
  arma::rowvec w(nrows);
  arma::vec y_hat(nrows);
  arma::vec residuals(nrows);

  for (int i =0; i< nrows; i++){
    w = kernel.row(i); //row vector
    // we have to convert w into a colvec; w.t()
    result = rcpp_coefficient(X, w.t(), Y);
    for (int c = 0; c < ncols; c++) {
      coef_result(i,c) = result(0, c);
    }
    // calculate predictions and residuals...
    y_hat.row(i) = X.row(i) * coef_result.row(i).t();
    residuals.row(i) = Y.row(i) - y_hat.row(i);
  }
  
  Rcpp::List res = Rcpp::List::create(
            Rcpp::Named("coef") = coef_result,
            Rcpp::Named("yfitted") = y_hat,
            Rcpp::Named("residuals") = residuals
            );
  
  return res; 
}

/*
 * Function for generating prediction. 
 * Not used. 
 */

//[[Rcpp::export]]
arma::mat rcpp_predicty(arma::mat coef, arma::mat X){
  int ncols = X.n_cols;
  int nrows = X.n_rows;
  // the output matrix for the results  
  arma::vec y_hat(nrows);
  arma::mat xrow(1, ncols);
  arma::mat cf(ncols, 1);
  double a;
  
  for (int i =0; i< nrows; i++){
    xrow = X.row(i); //row vector
    cf = coef.row(i); 
    // we have to convert w into a colvec; w.t()
    y_hat.row(i)= xrow * cf.t();
  }
  return y_hat;
}












