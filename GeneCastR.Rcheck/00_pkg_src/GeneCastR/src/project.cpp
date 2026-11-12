// Smoothing, projection, and the quantities read off them.
//
// Inside the observation window this is a smoother, not a filter: every
// survey in the deme, later ones included, informs the state at each step.
// Beyond the last survey the same kernel propagates the smoothed state
// forward, so the two halves of every curve join without a discontinuity.
//
// The forward quantity carried is the PREDICTIVE alpha^-_t -- after the
// transition, before the emission -- and the backward quantity beta_t
// includes the emission at t.  Pairing them counts each survey exactly once.
// Pairing the post-emission forward with the same backward would count the
// observed steps twice and make the posterior look sharper than it is.
//
// A deme carries one extra state, "not yet introduced", which emits a zero
// count with probability one and leaves at the t* hazard, depositing a single
// copy.  Because t* ~ U{1..T} has no support beyond the window, the mass
// still in that state at T can never become present: it stays in the
// denominator of every projected step.  That is a property of the model, not
// of this code, and it is why never-detected demes project low.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "GeneCastR.h"

static const double TINY = 1e-300;

// Smoothed lattice density and absent mass for one deme and one draw.
// `gam` is T x S row-major, `gp` is length T; both are normalised so that
// sum_k gam[t][k] + gp[t] == 1.
static void smooth_pt(const Banded& B, int N, int T, const DemeObs& o,
                      const std::vector<double>& cv, double tstar_p,
                      std::vector<double>& gam, std::vector<double>& gp) {
  const std::size_t S = (std::size_t)N + 1;

  // emissions by step
  std::vector<int> eidx((std::size_t)T + 1, -1);
  const int nobs = (int)o.step.size();
  std::vector<double> E((std::size_t)nobs * S);
  for (int j = 0; j < nobs; ++j) {
    if (o.step[j] >= 1 && o.step[j] <= T) {
      eidx[o.step[j]] = j;
      emission_fill(N, o.x[j], o.n[j], &E[(std::size_t)j * S]);
    }
  }

  std::vector<double> h;
  tstar_hazard(T, tstar_p, h);

  // ---- forward: predictive alpha^- and the absent mass ------------------
  std::vector<double> A((std::size_t)T * S), W((std::size_t)T);
  {
    std::vector<double> a(S, 0.0), nx(S, 0.0);
    double w = 1.0;
    for (int t = 1; t <= T; ++t) {
      const double inj = w * h[t - 1];
      w *= (1.0 - h[t - 1]);
      banded_tmv(B, a.data(), nx.data());
      for (std::size_t i = 0; i < S; ++i) nx[i] += inj * cv[i];
      a.swap(nx);
      double* Ar = &A[(std::size_t)(t - 1) * S];
      for (std::size_t i = 0; i < S; ++i) Ar[i] = a[i];
      W[t - 1] = w;
      const int j = eidx[t];
      if (j >= 0) {
        const double* e = &E[(std::size_t)j * S];
        double tot = 0.0;
        for (std::size_t i = 0; i < S; ++i) { a[i] *= e[i]; tot += a[i]; }
        if (o.x[j] != 0) w = 0.0;         // a positive count rules out absent
        tot += w;
        if (tot < TINY) tot = TINY;
        const double inv = 1.0 / tot;
        for (std::size_t i = 0; i < S; ++i) a[i] *= inv;
        w *= inv;
      }
    }
  }

  // ---- backward, combining into gamma as it goes ------------------------
  gam.assign((std::size_t)T * S, 0.0);
  gp.assign((std::size_t)T, 0.0);
  {
    std::vector<double> b(S, 1.0), nb(S, 0.0);
    double bp = 1.0;
    for (int t = T; t >= 1; --t) {
      if (t != T) {
        double cb = 0.0;
        for (std::size_t i = 0; i < S; ++i) cb += cv[i] * b[i];
        bp = h[t] * cb + (1.0 - h[t]) * bp;
        banded_mv(B, b.data(), nb.data());
        b.swap(nb);
      }
      const int j = eidx[t];
      if (j >= 0) {
        const double* e = &E[(std::size_t)j * S];
        for (std::size_t i = 0; i < S; ++i) b[i] *= e[i];
        if (o.x[j] != 0) bp = 0.0;
      }
      double sc = bp;
      for (std::size_t i = 0; i < S; ++i) if (b[i] > sc) sc = b[i];
      if (sc < TINY) sc = TINY;
      const double inv = 1.0 / sc;
      for (std::size_t i = 0; i < S; ++i) b[i] *= inv;
      bp *= inv;

      const double* Ar = &A[(std::size_t)(t - 1) * S];
      double* gr = &gam[(std::size_t)(t - 1) * S];
      double tot = W[t - 1] * bp;
      for (std::size_t i = 0; i < S; ++i) { gr[i] = Ar[i] * b[i]; tot += gr[i]; }
      if (tot < TINY) tot = TINY;
      const double g_inv = 1.0 / tot;
      for (std::size_t i = 0; i < S; ++i) gr[i] *= g_inv;
      gp[t - 1] = W[t - 1] * bp * g_inv;
    }
  }
}

//' Smoothed and projected quantities, averaged over posterior draws
//'
//' @param Z Matrix of log effective sizes, draws by demes.
//' @param U,NU Vectors of `u` and `nu`, one per draw.
//' @param M Lattice steps per unit time.
//' @param T Steps in the observation window.
//' @param H Steps to project beyond it.
//' @param demes Per-deme observations, as [GCData]`$deme_data()` returns.
//' @param thresholds Frequencies at which exceedance is reported.
//' @param nbins Bins for the frequency density.
//' @param founding `0` for a single copy at `t*`, `1` for a Poisson draw.
//' @param tstar_p Per-step founding hazard; non-positive keeps the uniform
//'   prior. A positive value lets arrival postdate the record, and the
//'   un-founded mass is carried into the projection at the same hazard.
//' @return A list with `pres` (the posterior mean presence probability),
//'   `pres_draws` (the per-draw trajectories the quantiles come from),
//'   `fmean`, `dens` (the mixture density over draws, excluding the atom at
//'   zero) and `exceed`.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List gc_project_pt(Rcpp::NumericMatrix Z, Rcpp::NumericVector U,
                         Rcpp::NumericVector NU, int M, int T, int H,
                         Rcpp::List demes, Rcpp::NumericVector thresholds,
                         int nbins, int founding = 0,
                         double tstar_p = -1.0) {
  const int D = demes.size();
  const int R = Z.nrow();
  const int TH = T + H;
  const int NT = thresholds.size();
  if (Z.ncol() != D) Rcpp::stop("Z must have one column per deme");
  if (U.size() != R || NU.size() != R)
    Rcpp::stop("U and NU must have one entry per row of Z");
  if (T < 1 || H < 0 || nbins < 1) Rcpp::stop("bad T, H or nbins");

  std::vector<DemeObs> obs((std::size_t)D);
  for (int d = 0; d < D; ++d) {
    Rcpp::List e = demes[d];
    Rcpp::IntegerVector st = e["step"], xx = e["x"], nn = e["n"];
    obs[d].step.assign(st.begin(), st.end());
    obs[d].x.assign(xx.begin(), xx.end());
    obs[d].n.assign(nn.begin(), nn.end());
    obs[d].cut = 0;
  }

  Rcpp::NumericMatrix pres(D, TH), fmean(D, TH);
  Rcpp::NumericVector pres_draws((std::size_t)R * D * TH);
  Rcpp::NumericVector dens((std::size_t)D * TH * nbins);
  Rcpp::NumericVector exceed((std::size_t)D * TH * (NT > 0 ? NT : 1));

  std::vector<double> gam, gp, p, np;

  for (int r = 0; r < R; ++r) {
    const double alpha = std::exp(U[r] / (double)M) - 1.0;
    for (int d = 0; d < D; ++d) {
      const double z = Z(r, d);
      double nv = (double)M * std::exp(z);
      int N = (int)std::nearbyint(nv);
      if (N < 2) N = 2;
      const std::size_t S = (std::size_t)N + 1;

      Banded B;
      const double lam = NU[r] * std::exp(z);
      banded_build(N, alpha, lam, B);
      std::vector<double> cv;
      founding_vec(N, lam, founding, cv);
      smooth_pt(B, N, T, obs[d], cv, tstar_p, gam, gp);

      // in-window
      for (int t = 0; t < T; ++t) {
        const double* g = &gam[(std::size_t)t * S];
        const double pr = 1.0 - gp[t] - g[0];
        pres(d, t) += pr;
        pres_draws[(std::size_t)r * D * TH + (std::size_t)d * TH + t] = pr;
        double fm = 0.0;
        for (std::size_t k = 1; k < S; ++k) {
          const double f = (double)k / (double)N;
          fm += f * g[k];
          int bin = (int)(f * (double)nbins);
          if (bin >= nbins) bin = nbins - 1;
          dens[((std::size_t)d * TH + t) * nbins + bin] += g[k];
        }
        fmean(d, t) += fm;
        for (int j = 0; j < NT; ++j) {
          const double thr = thresholds[j] * (double)N;
          double s = 0.0;
          for (std::size_t k = 1; k < S; ++k) if ((double)k > thr) s += g[k];
          exceed[((std::size_t)d * TH + t) * NT + j] += s;
        }
      }

      // projection: start from the SMOOTHED state at T, not a filtered one,
      // so the projected curve continues the in-window one without a step.
      //
      // Under a geometric t* some mass is still un-founded at T.  It is
      // carried forward at the same hazard and may be founded during the
      // horizon, which is the whole point of letting arrival postdate the
      // record: a deme that has not been reached yet can still be reached.
      p.assign(&gam[(std::size_t)(T - 1) * S], &gam[(std::size_t)T * S]);
      np.assign(S, 0.0);
      double w_un = gp[T - 1];
      const double hz = (tstar_p > 0.0) ? tstar_p : 0.0;
      for (int hh = 0; hh < H; ++hh) {
        banded_tmv(B, p.data(), np.data());
        const double inj = w_un * hz;
        w_un -= inj;
        if (inj > 0.0)
          for (std::size_t kk = 0; kk < S; ++kk) np[kk] += inj * cv[kk];
        p.swap(np);
        const int t = T + hh;
        double pr = 0.0, fm = 0.0;
        for (std::size_t k = 1; k < S; ++k) {
          const double f = (double)k / (double)N;
          pr += p[k];
          fm += f * p[k];
          int bin = (int)(f * (double)nbins);
          if (bin >= nbins) bin = nbins - 1;
          dens[((std::size_t)d * TH + t) * nbins + bin] += p[k];
        }
        pres(d, t) += pr;
        fmean(d, t) += fm;
        pres_draws[(std::size_t)r * D * TH + (std::size_t)d * TH + t] = pr;
        for (int j = 0; j < NT; ++j) {
          const double thr = thresholds[j] * (double)N;
          double s = 0.0;
          for (std::size_t k = 1; k < S; ++k) if ((double)k > thr) s += p[k];
          exceed[((std::size_t)d * TH + t) * NT + j] += s;
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }

  const double inv = 1.0 / (double)R;
  for (int i = 0; i < D * TH; ++i) { pres[i] *= inv; fmean[i] *= inv; }
  for (R_xlen_t i = 0; i < dens.size(); ++i) dens[i] *= inv;
  for (R_xlen_t i = 0; i < exceed.size(); ++i) exceed[i] *= inv;

  pres_draws.attr("dim") = Rcpp::IntegerVector::create(TH, D, R);
  dens.attr("dim")   = Rcpp::IntegerVector::create(nbins, TH, D);
  exceed.attr("dim") = Rcpp::IntegerVector::create(NT > 0 ? NT : 1, TH, D);

  return Rcpp::List::create(
    Rcpp::_["pres"] = pres, Rcpp::_["fmean"] = fmean,
    Rcpp::_["pres_draws"] = pres_draws,
    Rcpp::_["dens"] = dens, Rcpp::_["exceed"] = exceed);
}
