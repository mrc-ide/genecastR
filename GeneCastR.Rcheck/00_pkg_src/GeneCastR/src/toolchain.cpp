// Toolchain smoke test.  Confirms that Rcpp compiles and reports whether this
// build has OpenMP -- which Apple clang does not supply by default, so a
// package built on macOS without a separate libomp reports FALSE here and
// silently runs the per-deme loop serially.  Phase 7 depends on this being
// TRUE; everything before it is correct either way.
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

//' Report the compiled toolchain
//'
//' @return A list with `openmp` (logical), `max_threads` (integer, 1 when
//'   OpenMP is absent) and `cxx_standard` (the `__cplusplus` macro).
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List gc_toolchain() {
  bool has_omp = false;
  int threads = 1;
#ifdef _OPENMP
  has_omp = true;
  threads = omp_get_max_threads();
#endif
  return Rcpp::List::create(
    Rcpp::_["openmp"]       = has_omp,
    Rcpp::_["max_threads"]  = threads,
    Rcpp::_["cxx_standard"] = (double) __cplusplus);
}
