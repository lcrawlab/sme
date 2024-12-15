#include <Rcpp.h>

// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
bool check_openmp() {
  bool is_enabled = false;
#ifdef _OPENMP
  is_enabled = true;
#endif
  return is_enabled;
}
