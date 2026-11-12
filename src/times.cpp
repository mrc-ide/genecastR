// Forward-filter backward-sample, and the introduction times read off the
// sampled paths.
//
// T_first and T_est are PATH functionals: no marginal of the smoother
// determines them.  T_est asks when the excursion carrying the first
// DETECTED copy began -- walk back from the first survey with a non-zero
// count until the deme is empty -- which depends on the joint law of the
// path, so actual trajectories are needed.  It is undefined in a deme that
// was never detected, there being no observed lineage to date, and it is
// not conditioned on survival: the lineage that was seen need not be the
// one still present at T.  T_first is the exception -- in this model the founding event
// deposits a copy at t* and nothing is present before it, so T_first is
// identically t*, which the likelihood already marginalises analytically.
// Sampling it anyway is the check that the two routes agree.
//
// The sampler carries one extra state, N+1, meaning "not yet introduced".
// It emits a zero count with probability one, leaves at the t* hazard
// depositing a single copy, and is read back as k = 0 when the times are
// computed, so T_est is automatically at or after t*.  That row is the only
// time-varying part of the kernel and is applied analytically: building a
// kernel per step would be a (T, S, S) array, tens of gigabytes at the
// lattice sizes these chains visit.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "GeneCastR.h"

static const double TINY = 1e-300;

//' Introduction and establishment times, from sampled paths
//'
//' @param Z Matrix of log effective sizes, draws by demes.
//' @param U,NU Vectors of `u` and `nu`, one per draw.
//' @param M Lattice steps per unit time.
//' @param T Steps in the observation window.
//' @param demes Per-deme observations, as [GCData]`$deme_data()` returns.
//' @param n_paths Paths drawn per deme per draw.
//' @param joint Also return the joint distribution of the two times.
//' @param founding `0` for a single copy at `t*`, `1` for a Poisson draw.
//' @param tstar_p Per-step founding hazard; non-positive keeps the uniform
//'   prior. Under a geometric `t*` a path may still be un-founded at `T`,
//'   in which case it contributes to neither time and the rows of
//'   `pmf_tfirst` sum to less than one.
//' @return A list with `pmf_tfirst` and `pmf_test` (each `D x T`) and
//'   `p_alive`, the probability the variant is still present at `T`.
//'   `pmf_test` sums to one in a deme that was detected and to zero in one
//'   that was not, there being no observed lineage to date. `p_alive` is
//'   reported alongside but no longer conditions it.
//'   With `joint = TRUE`, a `T x T x D` array `joint` holding
//'   `p(T_first = i, T_est = j)`. Every path carries both times, so the
//'   pairing costs nothing beyond the storage; it is supported on
//'   `j >= i` only, since the surviving excursion cannot start before the
//'   variant first appears, and its diagonal is the probability that the
//'   lineage present at `T` is the one that originally arrived.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List gc_times_pt(Rcpp::NumericMatrix Z, Rcpp::NumericVector U,
                       Rcpp::NumericVector NU, int M, int T,
                       Rcpp::List demes, int n_paths,
                       bool joint = false, int founding = 0,
                       double tstar_p = -1.0,
                       Rcpp::Nullable<Rcpp::NumericVector> tstar_haz
                         = R_NilValue) {
  const int D = demes.size();
  const int R = Z.nrow();
  if (Z.ncol() != D) Rcpp::stop("Z must have one column per deme");
  if (U.size() != R || NU.size() != R)
    Rcpp::stop("U and NU must have one entry per row of Z");
  if (T < 1 || n_paths < 1) Rcpp::stop("bad T or n_paths");

  // T x D x R, per draw for the same reason the projection is
  const bool have_hz = tstar_haz.isNotNull();
  Rcpp::NumericVector HZ = have_hz ? Rcpp::NumericVector(tstar_haz.get())
                                   : Rcpp::NumericVector(0);
  if (have_hz) {
    Rcpp::IntegerVector dm = HZ.attr("dim");
    if (dm.size() != 3 || dm[0] < T || dm[1] != D || dm[2] != R)
      Rcpp::stop("tstar_haz must be T x D x nrow(Z)");
  }
  const std::size_t HZS = have_hz
    ? (std::size_t)Rcpp::as<Rcpp::IntegerVector>(HZ.attr("dim"))[0] : 0;

  std::vector<DemeObs> obs((std::size_t)D);
  for (int d = 0; d < D; ++d) {
    Rcpp::List e = demes[d];
    Rcpp::IntegerVector st = e["step"], xx = e["x"], nn = e["n"];
    obs[d].step.assign(st.begin(), st.end());
    obs[d].x.assign(xx.begin(), xx.end());
    obs[d].n.assign(nn.begin(), nn.end());
    // T_est needs the real cut -- the first step carrying a positive count --
    // because it asks when the excursion that produced THAT count began.
    Rcpp::IntegerVector cv = Rcpp::as<Rcpp::IntegerVector>(e["cut"]);
    obs[d].cut = (cv.size() == 0 || Rcpp::IntegerVector::is_na(cv[0]))
                 ? 0 : cv[0];
  }

  Rcpp::NumericMatrix pmf_first(D, T), pmf_est(D, T);
  Rcpp::NumericVector p_alive(D);
  Rcpp::NumericVector jnt(joint ? (R_xlen_t)T * T * D : 1);
  const double per = 1.0 / ((double)R * (double)n_paths);

  Rcpp::RNGScope scope;

  std::vector<double> A, h, a, nx, E, Pt, w;
  std::vector<int> kcur, tfirst, test;
  std::vector<char> alive, seen0;

  for (int r = 0; r < R; ++r) {
    const double alpha = std::exp(U[r] / (double)M) - 1.0;
    for (int d = 0; d < D; ++d) {
      const double z = Z(r, d);
      double nv = (double)M * std::exp(z);
      int N = (int)std::nearbyint(nv);
      if (N < 2) N = 2;
      const int SL = N + 1;        // lattice states
      const int S  = N + 2;        // plus "not yet introduced"

      Banded B;
      const double lamd = NU[r] * std::exp(z);
      banded_build(N, alpha, lamd, B);
      std::vector<double> cv;
      founding_vec(N, lamd, founding, cv);

      // dense column-major copy of the lattice kernel: the backward sample
      // needs a COLUMN of P per path, and a row-banded layout cannot give
      // one contiguously
      Pt.assign((std::size_t)SL * SL, 0.0);
      for (int i = 0; i < SL; ++i) {
        const double* row = &B.vals[B.off[i]];
        for (int jj = 0; jj < B.len[i]; ++jj)
          Pt[(std::size_t)(B.lo[i] + jj) * SL + i] = row[jj];
      }

      // emissions by step
      std::vector<int> eidx((std::size_t)T + 1, -1);
      const int nobs = (int)obs[d].step.size();
      E.assign((std::size_t)nobs * SL, 0.0);
      for (int j = 0; j < nobs; ++j)
        if (obs[d].step[j] >= 1 && obs[d].step[j] <= T) {
          eidx[obs[d].step[j]] = j;
          emission_fill(N, obs[d].x[j], obs[d].n[j], &E[(std::size_t)j * SL]);
        }

      if (have_hz) {
        const double* hp = &HZ[HZS * ((std::size_t)d + (std::size_t)D * r)];
        h.assign(hp, hp + T);
      } else tstar_hazard(T, tstar_p, h);

      // ---- forward filter ------------------------------------------------
      A.assign((std::size_t)(T + 1) * S, 0.0);
      a.assign(S, 0.0); nx.assign(S, 0.0);
      a[S - 1] = 1.0;
      for (int i = 0; i < S; ++i) A[i] = a[i];
      for (int t = 1; t <= T; ++t) {
        const double ht = h[t - 1], ab = a[S - 1];
        banded_tmv(B, a.data(), nx.data());        // lattice part only
        for (int i = 0; i < SL; ++i) nx[i] += ab * ht * cv[i];
        nx[S - 1] = ab * (1.0 - ht);
        const int j = eidx[t];
        if (j >= 0) {
          const double* e = &E[(std::size_t)j * SL];
          for (int i = 0; i < SL; ++i) nx[i] *= e[i];
          nx[S - 1] *= e[0];                       // absent emits like k = 0
        }
        double sm = 0.0;
        for (int i = 0; i < S; ++i) sm += nx[i];
        if (sm < TINY) sm = TINY;
        const double inv = 1.0 / sm;
        double* Ar = &A[(std::size_t)t * S];
        for (int i = 0; i < S; ++i) { nx[i] *= inv; Ar[i] = nx[i]; }
        a.swap(nx);
      }

      // ---- backward sample -----------------------------------------------
      kcur.assign(n_paths, 0);
      tfirst.assign(n_paths, -1);
      test.assign(n_paths, -1);
      alive.assign(n_paths, 0);
      seen0.assign(n_paths, 0);
      w.assign(S, 0.0);

      {
        const double* Ar = &A[(std::size_t)T * S];
        std::vector<double> cw(S);
        double acc = 0.0;
        for (int i = 0; i < S; ++i) { acc += Ar[i]; cw[i] = acc; }
        const int cut = obs[d].cut;            // first step with a positive count
        for (int i = 0; i < n_paths; ++i) {
          const double uu = R::unif_rand() * acc;
          int k = 0; while (k < S - 1 && cw[k] < uu) ++k;
          kcur[i] = k;
          const bool pos = (k >= 1 && k <= N);
          alive[i] = pos ? 1 : 0;                // reported, not a condition
          if (pos) tfirst[i] = T;
          // a deme that was never detected has no observed lineage to date
          seen0[i] = (cut <= 0) ? 1 : 0;
        }
      }

      for (int t = T - 1; t >= 0; --t) {
        const double ht = (t < T) ? h[t] : 0.0;
        const double* Ar = &A[(std::size_t)t * S];
        for (int i = 0; i < n_paths; ++i) {
          const int kn = kcur[i];
          double acc = 0.0;
          if (kn <= N) {
            const double* col = &Pt[(std::size_t)kn * SL];
            for (int j = 0; j < SL; ++j) { acc += Ar[j] * col[j]; w[j] = acc; }
            acc += Ar[S - 1] * ht * cv[kn];
            w[S - 1] = acc;
          } else {
            for (int j = 0; j < SL; ++j) w[j] = 0.0;
            acc = Ar[S - 1] * (1.0 - ht);
            w[S - 1] = acc;
          }
          if (acc < TINY) acc = TINY;
          const double uu = R::unif_rand() * acc;
          int k = 0; while (k < S - 1 && w[k] < uu) ++k;
          kcur[i] = k;
          const bool pos = (k >= 1 && k <= N);
          if (pos) tfirst[i] = t;
          // walking back from the first non-zero survey, the first empty
          // step we meet is the last one before it, so the excursion that
          // carried the detected copy began at t + 1
          if (!pos && !seen0[i] && t < obs[d].cut) {
            test[i] = t + 1;
            seen0[i] = 1;
          }
        }
        Rcpp::checkUserInterrupt();
      }

      for (int i = 0; i < n_paths; ++i) {
        if (tfirst[i] >= 1 && tfirst[i] <= T)
          pmf_first(d, tfirst[i] - 1) += per;
        if (alive[i]) p_alive[d] += per;
        if (test[i] >= 1 && test[i] <= T) {
          pmf_est(d, test[i] - 1) += per;
          if (joint && tfirst[i] >= 1)
            jnt[(R_xlen_t)d * T * T + (R_xlen_t)(test[i] - 1) * T
                + (tfirst[i] - 1)] += per;
        }
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["pmf_tfirst"] = pmf_first,
    Rcpp::_["pmf_test"]   = pmf_est,
    Rcpp::_["p_alive"]    = p_alive,
    Rcpp::_["joint"]      = R_NilValue);
  if (joint) {
    jnt.attr("dim") = Rcpp::IntegerVector::create(T, T, D);
    out["joint"] = jnt;
  }
  return out;
}
