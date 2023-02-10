#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Function to calculate the trace of S matrix
 *
 */
// [[Rcpp::export]]
arma::vec rcpp_Strace(const arma::mat X, const arma::mat kernel) {
  
  if (X.n_rows != kernel.n_cols) {
    throw std::runtime_error("kernel cols must equal rows of X");
  }
  
  if (X.n_rows != kernel.n_rows) {
    throw std::runtime_error("kernel rows must equal rows of X");
  }
  
  if (kernel.n_rows != kernel.n_cols) {
    throw std::runtime_error("kernel matrix must be square");
  }
  
  int nrows = X.n_rows;
  int ncols = kernel.n_cols;
  arma::colvec diagonal(ncols);
  arma::mat S(nrows, ncols);
  
  for (int i= 0; i < nrows; i++){
    
    for (int j = 0; j < ncols; j++){
      diagonal(j) = kernel(i, j);
    }
    
    S.row(i) = X.row(i) * inv(X.t()*diagmat(diagonal)*X) * X.t() * diagmat(diagonal);
  }
  
  return S.diag();
  
}


/*
 * simple function for standard deviation; just to make it more clear
 *
 */

// [[Rcpp::export]]
double armaSD(arma::colvec inVec)
{
  return arma::stddev(inVec);
}


/*
 * function for calculating for calculating AICc
 *
 */


// [[Rcpp::export]]
double rcpp_AICc(const arma::mat X, const arma::mat kernel, const arma::colvec residuals) {
  
  double aicc;
  double n = X.n_rows;
  double sigma_hat = armaSD(residuals); //checked - same as R's st(df$residuals)
  arma::vec St = rcpp_Strace(X, kernel); // checked - same as R's diag(P)

  double sumSt = arma::accu(St); // checked - same as R's sum(diag(P))
  
  // M_PI - this is the same as 'pi' //checked
  aicc = (2.0*n*std::log(sigma_hat)) + (n*std::log(2.0*M_PI)) + (n*((n+sumSt)/(n-2-sumSt)));
  
  return aicc;
  
}


