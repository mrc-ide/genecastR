#ifndef GENECASTR_H
#define GENECASTR_H

#include <vector>
#include <cstddef>

// A row-banded transition kernel.  Row i holds `len[i]` live entries starting
// at column `lo[i]`, packed contiguously at `vals[off[i]]`.  Everything
// outside the band is below 1e-18 of the row's maximum, and each row is
// renormalised so it stays exactly stochastic.  Row 0 -- the arrival row out
// of the empty state -- is kept dense, because its folded Poisson tail
// reaches the last cell.
struct Banded {
  int N;
  std::vector<int> lo, len;
  std::vector<std::size_t> off;
  std::vector<double> vals;
};

void banded_build(int N, double alpha, double lam, Banded& B);

// The law of the count at t*.
//
//   founding = 0  a single copy, exactly: the deme is founded by one
//                 immigrant and the flow process governs only what happens
//                 afterwards.
//   founding = 1  Poisson(lambda), the k = 0 cell included: t* is the step
//                 from which the arrival RATE is positive, so the deme
//                 joins the flow process then and its first actual arrival
//                 happens at or after t*.  This is the no-t* model with the
//                 deme's start time moved from step 1 to t*.
void founding_vec(int N, double lam, int founding, std::vector<double>& cv);
void poisson_row_fill(int N, double lam, double* row);
void emission_fill(int N, int x, int n, double* e);
double gc_dot8(const double* a, const double* b, std::size_t n);

// out = P b, over each row's live band.
inline void banded_mv(const Banded& B, const double* b, double* out) {
  const std::size_t S = (std::size_t)B.N + 1;
  const double* v = B.vals.data();
  for (std::size_t i = 0; i < S; ++i)
    out[i] = gc_dot8(v + B.off[i], b + B.lo[i], (std::size_t)B.len[i]);
}

// out = a P, the transpose product, as a scatter over the same bands.
inline void banded_tmv(const Banded& B, const double* a, double* out) {
  const std::size_t S = (std::size_t)B.N + 1;
  for (std::size_t j = 0; j < S; ++j) out[j] = 0.0;
  const double* v = B.vals.data();
  for (std::size_t i = 0; i < S; ++i) {
    const double ai = a[i];
    if (ai == 0.0) continue;
    const double* row = v + B.off[i];
    double* dst = out + B.lo[i];
    const int n = B.len[i];
    for (int j = 0; j < n; ++j) dst[j] += ai * row[j];
  }
}

// Lambda[t-1] = log p(X | t* = t, c = 1), for t = 1..T, with -Inf at and
// after `cut`.  `cut <= 0` means the deme was never positive, so t* is
// unconstrained.
// Lambda[t-1] = log p(X | t* = t) + log p(t* = t), for t = 1..T, and
// Lambda[T] = log p(X | t* > T) + log p(t* > T) -- the branch in which the
// variant has not reached the deme within the record at all.  `out` must
// therefore hold T + 1 doubles.
//
//   tstar_p <= 0   t* ~ Uniform{1..T}; the deme is founded inside the
//                  window with probability one and the last element is -Inf.
//   tstar_p  > 0   t* ~ Geometric(tstar_p) per step, so arrival carries a
//                  constant hazard and may fall beyond the record.
void lambda_pt_fill(int N, double alpha, double lam, int T,
                    const int* step, const int* x, const int* n,
                    int nobs, int cut, double* out, int founding = 0,
                    double tstar_p = -1.0);

double gc_logsumexp(const double* v, int m);

// One deme's observations, copied out of R once so an inner loop never
// touches an R object.
// Hazard of founding at each step, given it has not happened yet.
//
//   tstar_p <= 0   t* ~ Uniform{1..T}: the hazard rises to one at the last
//                  step, so the deme is certainly founded inside the window.
//   tstar_p  > 0   t* ~ Geometric: a constant hazard, so (1-p)^T of the mass
//                  is still un-founded at T and arrival may postdate the
//                  record entirely.
inline void tstar_hazard(int T, double tstar_p, std::vector<double>& h) {
  h.assign((std::size_t)T, 0.0);
  if (tstar_p > 0.0) {
    for (int t = 0; t < T; ++t) h[t] = tstar_p;
    return;
  }
  double rem = 1.0;
  const double pri = 1.0 / (double)T;
  for (int t = 0; t < T; ++t) {
    h[t] = (rem > 1e-300) ? pri / rem : 0.0;
    rem -= pri;
  }
}

struct DemeObs {
  std::vector<int> step, x, n;
  int cut;
};

#endif
