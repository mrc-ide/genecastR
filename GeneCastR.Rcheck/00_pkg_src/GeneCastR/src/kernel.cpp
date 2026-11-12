// Wright-Fisher lattice kernels, the t*-marginalised backward pass, and the
// quantities read off it: the per-deme log-likelihood and the t* posterior.
//
// The backward recursion is
//
//     beta_t(k) = p(X_{obs >= t} | k_t = k),   beta_t = e_t * (P beta_{t+1}),
//
// so an introduction of one copy at t* has likelihood exactly beta_{t*}(1),
// provided every observation before t* is zero -- which it must be, because
// the frequency is identically zero then.  One backward pass therefore yields
// the whole t* likelihood surface.
//
// Only column 1 of that surface is ever used.  The reference implementation
// stores the full (T, N+1) table and takes its logarithm; here a single
// scalar per step is kept, which removes T*(N+1) logarithms per deme per
// call.  The per-step rescaling and its accumulated logarithm are reproduced
// exactly, floor at 1e-300 included, so the two agree to rounding.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
#include "GeneCastR.h"

static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double TINY = 1e-300;

// A dot product written as eight independent accumulators.
//
// The obvious single-accumulator loop does not vectorise: floating-point
// addition is not associative, so the compiler may not reorder the reduction
// without -ffast-math, and every step then waits on the previous add's
// latency.  On this lattice that is the whole cost of the backward pass --
// the naive form ran the Ugandan panel seven times slower than the reference
// NumPy implementation, which gets reassociation for free inside BLAS.
// Splitting the chain restores it with no change to the compiler flags.
//
// This does reassociate the sum, so results differ from any other summation
// order in the last few bits.  That is already true of BLAS, and is far
// below the 1e-10 at which the parity tests compare.
double gc_dot8(const double* __restrict a,
                          const double* __restrict b, std::size_t n) {
  double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0;
  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    s0 += a[i    ] * b[i    ];
    s1 += a[i + 1] * b[i + 1];
    s2 += a[i + 2] * b[i + 2];
    s3 += a[i + 3] * b[i + 3];
    s4 += a[i + 4] * b[i + 4];
    s5 += a[i + 5] * b[i + 5];
    s6 += a[i + 6] * b[i + 6];
    s7 += a[i + 7] * b[i + 7];
  }
  for (; i < n; ++i) s0 += a[i] * b[i];
  return ((s0 + s1) + (s2 + s3)) + ((s4 + s5) + (s6 + s7));
}

// Genic selection: psi(f, alpha) = f (1 + alpha) / (1 + alpha f).
static inline double psi(double f, double alpha) {
  return f * (1.0 + alpha) / (1.0 + alpha * f);
}

// ----------------------------------------------------------- banded form --

// A binomial row is numerically zero away from its mode: at N = 180 the
// standard deviation is about 6 copies, so fewer than half the entries carry
// anything a double can represent relative to the row maximum.  Storing and
// multiplying only the live band halves both the arithmetic and the memory
// traffic of the backward pass, and halves the number of exponentials in the
// kernel build.
//
// The cut is placed where the log pmf has fallen 41.4 below its mode, a ratio
// of 1e-18, and each row is renormalised afterwards so it remains exactly
// stochastic -- mass conservation is what keeps a 112-step backward pass
// stable, and is worth more than the 1e-18 that was dropped.

void banded_build(int N, double alpha, double lam, Banded& B) {
  const int S = N + 1;
  B.N = N;
  B.lo.assign(S, 0);
  B.len.assign(S, 0);
  B.off.assign(S, 0);
  B.vals.clear();
  B.vals.reserve((std::size_t)S * 48);

  std::vector<double> lch(S), scratch(S);
  const double lgN = std::lgamma((double)N + 1.0);
  for (int j = 0; j < S; ++j)
    lch[j] = lgN - std::lgamma((double)j + 1.0)
                 - std::lgamma((double)(N - j) + 1.0);

  // row 0 is the Poisson arrival row, kept dense: its tail fold reaches N
  {
    std::vector<double> r(S);
    poisson_row_fill(N, lam, &r[0]);
    B.lo[0] = 0; B.len[0] = S; B.off[0] = 0;
    B.vals.insert(B.vals.end(), r.begin(), r.end());
  }

  for (int k = 1; k < S; ++k) {
    double p = psi((double)k / (double)N, alpha);
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    B.off[k] = B.vals.size();
    if (p >= 1.0) { B.lo[k] = N; B.len[k] = 1; B.vals.push_back(1.0); continue; }
    if (p <= 0.0) { B.lo[k] = 0; B.len[k] = 1; B.vals.push_back(1.0); continue; }

    // One exponential at the mode, then the ratio recursion outwards:
    //
    //     pmf(j+1) / pmf(j) = ((N - j) / (j + 1)) * p / (1 - p),
    //
    // so each further entry is a single multiply.  Evaluating the whole band
    // through exp() instead costs an exponential per entry and dominated the
    // kernel build.  Walking outwards from the mode the values only decrease,
    // so the relative error does not grow, and the row is renormalised at the
    // end in any case.
    const double lp = std::log(p), l1p = std::log1p(-p);
    const double q = p / (1.0 - p);
    int m = (int)((double)(N + 1) * p);
    if (m < 0) m = 0;
    if (m > N) m = N;
    const double vm = std::exp(lch[m] + (double)m * lp + (double)(N - m) * l1p);
    const double thr = vm * 1e-18;
    scratch[m] = vm;

    int hi = m;
    double v = vm;
    for (int j = m; j < N; ++j) {
      v *= ((double)(N - j) / (double)(j + 1)) * q;
      if (!(v > thr) || !std::isfinite(v)) break;
      scratch[j + 1] = v;
      hi = j + 1;
    }
    int lo = m;
    v = vm;
    for (int j = m; j > 0; --j) {
      v *= ((double)j / (double)(N - j + 1)) / q;
      if (!(v > thr) || !std::isfinite(v)) break;
      scratch[j - 1] = v;
      lo = j - 1;
    }

    const int len = hi - lo + 1;
    B.lo[k] = lo; B.len[k] = len;
    double sum = 0.0;
    for (int j = lo; j <= hi; ++j) sum += scratch[j];
    const double inv = (sum > 0.0) ? 1.0 / sum : 1.0;
    const std::size_t at = B.vals.size();
    B.vals.resize(at + (std::size_t)len);
    double* dst = &B.vals[at];
    for (int j = 0; j < len; ++j) dst[j] = scratch[lo + j] * inv;
  }
}

// ---------------------------------------------------------------- kernel --

// Row k of the Wright-Fisher kernel is Binomial(N, psi(k/N, alpha)).  Rows 0
// and N are absorbing: psi is exactly 0 and 1 there, and the general formula
// would need 0*log(0) = 0 to reach the same answer, so they are written out.
static void gc_kernel_fill(int N, double alpha, std::vector<double>& P) {
  const std::size_t S = (std::size_t)N + 1;
  P.assign(S * S, 0.0);

  std::vector<double> lch(S);
  const double lgN = std::lgamma((double)N + 1.0);
  for (int j = 0; j <= N; ++j) {
    lch[j] = lgN - std::lgamma((double)j + 1.0)
                 - std::lgamma((double)(N - j) + 1.0);
  }

  for (int k = 0; k <= N; ++k) {
    double p = psi((double)k / (double)N, alpha);
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    double* row = &P[(std::size_t)k * S];
    if (p <= 0.0) { row[0] = 1.0; continue; }
    if (p >= 1.0) { row[N] = 1.0; continue; }
    const double lp = std::log(p), l1p = std::log1p(-p);
    for (int j = 0; j <= N; ++j) {
      row[j] = std::exp(lch[j] + (double)j * lp + (double)(N - j) * l1p);
    }
  }
}

// Poisson arrivals out of the empty state, truncated at N with the tail
// folded onto the last cell so the row sums to exactly one.  Leaking even
// 1e-12 here biases the likelihood downwards over T steps, and it bites
// hardest in the small demes where the tail is heaviest relative to N.
void poisson_row_fill(int N, double lam, double* row) {
  if (lam <= 0.0) {
    for (int k = 0; k <= N; ++k) row[k] = 0.0;
    row[0] = 1.0;
    return;
  }
  const double loglam = std::log(lam);
  double s = 0.0;
  for (int k = 0; k <= N; ++k) {
    const double lt = -lam + (k == 0 ? 0.0 : (double)k * loglam)
                      - std::lgamma((double)k + 1.0);
    row[k] = std::exp(lt);
    s += row[k];
  }
  const double rem = 1.0 - s;
  row[N] += (rem > 0.0 ? rem : 0.0);
}

void founding_vec(int N, double lam, int founding, std::vector<double>& cv) {
  cv.assign((std::size_t)N + 1, 0.0);
  if (founding == 0) {
    cv[1] = 1.0;                       // exactly one copy
    return;
  }
  poisson_row_fill(N, lam, &cv[0]);    // the law that governs re-entry
  if (founding == 1) return;           // t* = onset of the rate; 0 allowed

  // founding == 2: t* is the first NON-ZERO draw, so the founding size is
  // Poisson(lambda) conditioned on being an arrival at all.
  const double p0 = cv[0];
  cv[0] = 0.0;
  const double rest = 1.0 - p0;
  if (rest > 1e-12) {
    const double inv = 1.0 / rest;
    for (int i = 1; i <= N; ++i) cv[i] *= inv;
  } else {
    // lambda so small that the conditional law is numerically a point mass
    // at one copy, which is also its exact limit as lambda -> 0
    for (int i = 1; i <= N; ++i) cv[i] = 0.0;
    cv[1] = 1.0;
  }
}

// Binomial(n, k/N) evaluated at x, over the whole lattice.
void emission_fill(int N, int x, int n, double* e) {
  const double lg = std::lgamma((double)n + 1.0)
                  - std::lgamma((double)x + 1.0)
                  - std::lgamma((double)(n - x) + 1.0);
  for (int k = 0; k <= N; ++k) {
    const double f = (double)k / (double)N;
    double a, b;
    if (x == 0)          a = 0.0;
    else if (f > 0.0)    a = (double)x * std::log(f);
    else                 a = NEG_INF;
    if (n - x == 0)      b = 0.0;
    else if (f < 1.0)    b = (double)(n - x) * std::log1p(-f);
    else                 b = NEG_INF;
    e[k] = std::exp(lg + a + b);
  }
}

// ----------------------------------------------------------- the Lambda --

// Lambda[t - 1] = log p(X | t* = t, c = 1), for t = 1..T, with -Inf at and
// after `cut` (the first step carrying a non-zero count): the variant cannot
// arrive after it has already been seen.  `cut <= 0` means no constraint.
void lambda_pt_fill(int N, double alpha, double lam, int T,
                    const int* step, const int* x, const int* n,
                    int nobs, int cut, double* out, int founding,
                    double tstar_p) {
  const std::size_t S = (std::size_t)N + 1;

  Banded B;
  banded_build(N, alpha, lam, B);
  std::vector<double> cv;
  founding_vec(N, lam, founding, cv);

  // emissions, indexed by step, built once
  std::vector<int> eidx((std::size_t)T + 1, -1);
  std::vector<double> E((std::size_t)nobs * S);
  for (int j = 0; j < nobs; ++j) {
    if (step[j] >= 1 && step[j] <= T) {
      eidx[step[j]] = j;
      emission_fill(N, x[j], n[j], &E[(std::size_t)j * S]);
    }
  }

  std::vector<double> beta(S, 1.0), tmp(S, 0.0);
  double acc = 0.0;

  for (int t = T; t >= 1; --t) {
    if (t != T) {
      // tmp = P beta, over the live band of each row only
      const double* vp = &B.vals[0];
      const double* bp = &beta[0];
      for (std::size_t i = 0; i < S; ++i)
        tmp[i] = gc_dot8(vp + B.off[i], bp + B.lo[i], (std::size_t)B.len[i]);
      beta.swap(tmp);
    }
    const int j = eidx[t];
    if (j >= 0) {
      const double* e = &E[(std::size_t)j * S];
      for (std::size_t i = 0; i < S; ++i) beta[i] *= e[i];
    }
    double mx = 0.0;
    for (std::size_t i = 0; i < S; ++i) if (beta[i] > mx) mx = beta[i];
    if (mx < TINY) mx = TINY;
    const double inv = 1.0 / mx;
    for (std::size_t i = 0; i < S; ++i) beta[i] *= inv;
    acc += std::log(mx);
    // integrate the founding law against the backward message
    double w = 0.0;
    for (std::size_t i = 0; i < S; ++i) w += cv[i] * beta[i];
    out[t - 1] = std::log(w < TINY ? TINY : w) + acc;
  }

  // the prior over t*
  if (tstar_p > 0.0) {
    // Geometric: p(t* = t) = (1-p)^(t-1) p, and the remaining (1-p)^T is the
    // branch where the variant never reached the deme within the record.
    const double l1p = std::log1p(-tstar_p), lp = std::log(tstar_p);
    for (int t = 1; t <= T; ++t) out[t - 1] += (double)(t - 1) * l1p + lp;
    // p(X | never founded) is one when every count is zero and zero
    // otherwise: an empty deme cannot produce a positive observation.
    bool all_zero = true;
    for (int j = 0; j < nobs; ++j) if (x[j] != 0) { all_zero = false; break; }
    out[T] = all_zero ? (double)T * l1p : NEG_INF;
  } else {
    const double lu = -std::log((double)T);
    for (int t = 1; t <= T; ++t) out[t - 1] += lu;
    out[T] = NEG_INF;                  // uniform: arrival inside the window
  }

  if (cut > 0 && cut < T) {
    for (int t = cut + 1; t <= T; ++t) out[t - 1] = NEG_INF;
  }
}

double gc_logsumexp(const double* v, int m) {
  double mx = NEG_INF;
  for (int i = 0; i < m; ++i) if (v[i] > mx) mx = v[i];
  if (!std::isfinite(mx)) return mx;
  double s = 0.0;
  for (int i = 0; i < m; ++i) s += std::exp(v[i] - mx);
  return mx + std::log(s);
}

// ------------------------------------------------------------- exported --

//' Wright-Fisher transition kernel on the lattice
//'
//' @param N Lattice size; the state space is `0..N`.
//' @param alpha Per-step selection coefficient.
//' @return An `(N+1) x (N+1)` row-stochastic matrix. Rows 0 and N are
//'   absorbing.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix gc_kernel(int N, double alpha) {
  if (N < 1) Rcpp::stop("N must be at least 1");
  if (alpha <= -1.0) Rcpp::stop("alpha must exceed -1");
  std::vector<double> P;
  gc_kernel_fill(N, alpha, P);
  Rcpp::NumericMatrix out(N + 1, N + 1);
  for (int i = 0; i <= N; ++i)
    for (int j = 0; j <= N; ++j)
      out(i, j) = P[(std::size_t)i * (N + 1) + j];
  return out;
}

//' Poisson arrival row for the empty state
//'
//' @param N Lattice size.
//' @param lambda Per-step Poisson mean.
//' @return A length `N+1` probability vector summing to exactly 1, the tail
//'   beyond `N` folded onto the last cell.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gc_poisson_row(int N, double lambda) {
  if (N < 1) Rcpp::stop("N must be at least 1");
  if (lambda < 0.0) Rcpp::stop("lambda must be non-negative");
  Rcpp::NumericVector out(N + 1);
  poisson_row_fill(N, lambda, &out[0]);
  return out;
}

//' Binomial emission over the lattice
//'
//' @param N Lattice size.
//' @param x,n Observed count and sample size.
//' @return `p(x | k, n)` for `k = 0..N`.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gc_emission(int N, int x, int n) {
  if (N < 1) Rcpp::stop("N must be at least 1");
  if (x < 0 || n < x) Rcpp::stop("need 0 <= x <= n");
  Rcpp::NumericVector out(N + 1);
  emission_fill(N, x, n, &out[0]);
  return out;
}

// Pull one deme's observations out of the list element.
static void deme_fields(const Rcpp::List& e, Rcpp::IntegerVector& step,
                        Rcpp::IntegerVector& x, Rcpp::IntegerVector& n,
                        int& cut) {
  step = Rcpp::as<Rcpp::IntegerVector>(e["step"]);
  x    = Rcpp::as<Rcpp::IntegerVector>(e["x"]);
  n    = Rcpp::as<Rcpp::IntegerVector>(e["n"]);
  SEXP cs = e["cut"];
  Rcpp::IntegerVector cv = Rcpp::as<Rcpp::IntegerVector>(cs);
  cut = (cv.size() == 0 || Rcpp::IntegerVector::is_na(cv[0])) ? 0 : cv[0];
  if (step.size() != x.size() || step.size() != n.size())
    Rcpp::stop("step, x and n must have the same length");
}

//' Per-deme log-likelihood, introduction time marginalised
//'
//' @param Nlat Integer vector of lattice sizes, one per deme.
//' @param alpha Per-step selection coefficient, shared.
//' @param lambda Numeric vector of per-step Poisson means, one per deme.
//' @param T Number of lattice steps.
//' @param demes List of per-deme observations, as [GCData]`$deme_data()`
//'   returns: `step`, `x`, `n` and `cut`.
//' @param founding `0` for a single copy deposited at `t*`, `1` for a
//'   Poisson draw, which makes `t*` the step from which the arrival rate is
//'   positive rather than the arrival itself.
//' @param tstar_p Per-step founding hazard. Non-positive keeps the uniform
//'   prior on `t*`, under which the deme is founded inside the window with
//'   probability one. A positive value makes `t*` geometric, so arrival may
//'   fall beyond the record.
//' @return `log p(X_d | theta)` for each deme, with `t*` uniform on `1..T`
//'   and integrated out.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gc_loglik_pt(Rcpp::IntegerVector Nlat, double alpha,
                                 Rcpp::NumericVector lambda, int T,
                                 Rcpp::List demes, int founding = 0,
                                 double tstar_p = -1.0) {
  const int D = demes.size();
  if (Nlat.size() != D || lambda.size() != D)
    Rcpp::stop("Nlat and lambda must have one entry per deme");
  if (T < 1) Rcpp::stop("T must be at least 1");

  Rcpp::NumericVector out(D);
  std::vector<double> lam((std::size_t)T + 1);

  for (int d = 0; d < D; ++d) {
    Rcpp::List e = demes[d];
    Rcpp::IntegerVector step, x, n;
    int cut;
    deme_fields(e, step, x, n, cut);
    lambda_pt_fill(Nlat[d], alpha, lambda[d], T, &step[0], &x[0], &n[0],
                   step.size(), cut, &lam[0], founding, tstar_p);
    out[d] = gc_logsumexp(&lam[0], T + 1);
  }
  return out;
}

//' Posterior over the introduction time, for one parameter draw
//'
//' Analytic: a normalisation of the same backward row the likelihood sums.
//' The introduction time is never instantiated, so this carries no
//' Monte-Carlo error.
//'
//' @inheritParams gc_loglik_pt
//' @return A `D x (T+1)` matrix; row `d` is `p(t* = t | X_d, theta)` for
//'   `t = 1..T` with a final column holding `p(t* > T | X_d, theta)`, the
//'   probability the variant had not reached the deme within the record.
//'   The row sums to one. Under the uniform prior that last column is zero.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix gc_tstar_pt(Rcpp::IntegerVector Nlat, double alpha,
                                Rcpp::NumericVector lambda, int T,
                                Rcpp::List demes, int founding = 0,
                                double tstar_p = -1.0) {
  const int D = demes.size();
  if (Nlat.size() != D || lambda.size() != D)
    Rcpp::stop("Nlat and lambda must have one entry per deme");
  if (T < 1) Rcpp::stop("T must be at least 1");

  Rcpp::NumericMatrix out(D, T + 1);
  std::vector<double> lam((std::size_t)T + 1);

  for (int d = 0; d < D; ++d) {
    Rcpp::List e = demes[d];
    Rcpp::IntegerVector step, x, n;
    int cut;
    deme_fields(e, step, x, n, cut);
    lambda_pt_fill(Nlat[d], alpha, lambda[d], T, &step[0], &x[0], &n[0],
                   step.size(), cut, &lam[0], founding, tstar_p);
    const double lse = gc_logsumexp(&lam[0], T + 1);
    for (int t = 0; t <= T; ++t) out(d, t) = std::exp(lam[t] - lse);
  }
  return out;
}

//' Banded kernel statistics
//'
//' Diagnostic: how much of the dense kernel is numerically live.
//'
//' @param N Lattice size.
//' @param alpha Per-step selection coefficient.
//' @param lambda Per-step Poisson mean for the empty state.
//' @return A list with the mean and maximum band width, the dense width, and
//'   the fraction of entries retained.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List gc_band_stats(int N, double alpha, double lambda) {
  Banded B;
  banded_build(N, alpha, lambda, B);
  double tot = 0.0; int mx = 0;
  for (int i = 0; i <= N; ++i) { tot += B.len[i]; if (B.len[i] > mx) mx = B.len[i]; }
  return Rcpp::List::create(
    Rcpp::_["mean_width"] = tot / (double)(N + 1),
    Rcpp::_["max_width"]  = mx,
    Rcpp::_["dense"]      = N + 1,
    Rcpp::_["kept"]       = tot / ((double)(N + 1) * (double)(N + 1)));
}
