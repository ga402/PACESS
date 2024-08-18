#include <Rcpp.h>
using namespace Rcpp;

// This is a function for generating distance matrix for 3D data.

// [[Rcpp::export]]
NumericMatrix rcpp_distance3d(NumericMatrix mat){
  
  if (mat.ncol() != 3) {
    throw std::runtime_error("Incompatible number of dimensions");
  }

  int nrows = mat.nrow();
  int ncols = mat.ncol();
  
  NumericMatrix dmat(nrows, nrows);
  
  for (int r1 = 0; r1 < nrows; r1++) {
    for (int r2 = 0; r2 < nrows; r2++) {
      double total = 0;
      for (int c = 0; c < ncols; c++) {
        total += std::pow(mat(r1, c) - mat(r2, c), 2);
      }
      dmat(r1, r2) = std::sqrt(total);
    }
  }
  
  return dmat;
}



