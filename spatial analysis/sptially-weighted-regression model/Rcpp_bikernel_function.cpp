#include <Rcpp.h>
using namespace Rcpp;

/*
 * Function to calculate the bikernel function
 * 
 */

// [[Rcpp::export]]
NumericMatrix rcpp_bikernel(NumericMatrix dmat, double b){

  if (dmat.ncol() != dmat.nrow()) {
    throw std::runtime_error("dmat should be a square matrix");
  }
  
  int nrows = dmat.nrow();
  
  NumericMatrix bikernel_mat(nrows, nrows);
  
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) {
      if (dmat(i,j) < b){
        bikernel_mat(i, j) = std::pow(1-std::pow((dmat(i,j)/b), 2),2);
      } else{
        bikernel_mat(i,j) = 0;
      }
    }
  }
  return bikernel_mat;
}



