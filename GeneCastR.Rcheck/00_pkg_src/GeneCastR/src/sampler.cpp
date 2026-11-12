// Blocked Metropolis-within-Gibbs for the hierarchical model.
//
// The likelihood is exact, so this is ordinary Metropolis-within-Gibbs and
// not a pseudo-marginal scheme: the acceptance ratio uses the true value and
// the chain targets the posterior without any auxiliary-variable argument.
//
// Four blocks per sweep:
//
//   1. every z_d at once.  They are conditionally independent given the
//      shared parameters, so one vectorised pass proposes all D and each is
//      accepted or rejected against its own carried log-likelihood.
//   2. u = log(1 + A), random walk, Gaussian prior.
//   3. log nu, random walk, half-normal prior with the Jacobian.
//   4. mu and sigma^2 from their exact full conditionals -- no likelihood.
//
// The N_eff and u bounds are rejections applied BEFORE the likelihood is
// evaluated.  That is not only a saving: it is what stops an absurd proposal
// allocating a vast lattice, and it supplies the negative feedback that keeps
// step-size adaptation stable.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "GeneCastR.h"

static inline double loglik_one(const DemeObs& o, int N, double alpha,
                                double lam, int T, std::vector<double>& buf,
                                int founding, double tstar_p) {
  lambda_pt_fill(N, alpha, lam, T, o.step.data(), o.x.data(), o.n.data(),
                 (int)o.step.size(), o.cut, buf.data(), founding, tstar_p);
  return gc_logsumexp(buf.data(), T + 1);
}

static inline int nlat_of(double z, int M) {
  double v = (double)M * std::exp(z);
  // R's round(), which is half-to-even; std::nearbyint with the default
  // rounding mode does the same, and matches gc_nlat() exactly.
  int n = (int)std::nearbyint(v);
  return n < 2 ? 2 : n;
}

//' Blocked Metropolis-within-Gibbs for the Poisson model with t*
//'
//' @param demes Per-deme observations, as [GCData]`$deme_data()` returns.
//' @param M Lattice steps per unit time.
//' @param T Number of lattice steps.
//' @param n_iter Sweeps to run, warm-up included.
//' @param warm Warm-up sweeps; step sizes adapt only within these.
//' @param z0 Starting log effective sizes, one per deme.
//' @param u0,nu0 Starting values for the shared parameters.
//' @param step_z,step_u,step_nu Starting random-walk scales.
//' @param prior Named list of prior settings.
//' @param target Acceptance rate the adaptation aims at.
//' @param founding `0` for a single copy at `t*`, `1` for a Poisson draw.
//' @return A list with the draws, the acceptance rates, the final step sizes
//'   and the count of proposals rejected by the bounds before any likelihood
//'   was evaluated.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List gc_mwg_run(Rcpp::List demes, int M, int T, int n_iter, int warm,
                      Rcpp::NumericVector z0, double u0, double nu0,
                      double step_z, double step_u, double step_nu,
                      Rcpp::List prior, double target, int founding = 0,
                      double tstar_p = -1.0, bool shared = false,
                      int u_mode = 0) {

  const int D = demes.size();
  if (z0.size() != D) Rcpp::stop("z0 must have one entry per deme");
  if (n_iter < 1) Rcpp::stop("n_iter must be at least 1");
  if (warm < 0 || warm >= n_iter) Rcpp::stop("need 0 <= warm < n_iter");

  // copy the data out of R once
  std::vector<DemeObs> obs((std::size_t)D);
  for (int d = 0; d < D; ++d) {
    Rcpp::List e = demes[d];
    Rcpp::IntegerVector st = e["step"], xx = e["x"], nn = e["n"];
    Rcpp::IntegerVector cv = Rcpp::as<Rcpp::IntegerVector>(e["cut"]);
    obs[d].step.assign(st.begin(), st.end());
    obs[d].x.assign(xx.begin(), xx.end());
    obs[d].n.assign(nn.begin(), nn.end());
    obs[d].cut = (cv.size() == 0 || Rcpp::IntegerVector::is_na(cv[0]))
                 ? 0 : cv[0];
  }

  const double tau    = prior["tau"];
  const double u_lo   = prior["u_lo"],   u_hi   = prior["u_hi"];
  const double nu_s   = prior["nu_s"];
  const double lnu_lo = prior["lnu_lo"], lnu_hi = prior["lnu_hi"];
  const double mu0    = prior["mu0"],    s0     = prior["s0"];
  const double ig_a   = prior["ig_a"],   ig_b   = prior["ig_b"];
  const double ne_lo  = prior["ne_lo"],  ne_hi  = prior["ne_hi"];
  // theta: the per-YEAR founding hazard, when t* is geometric.  Its prior
  // is Exponential(th_rate); with th_rate set to the window length in years
  // that makes P(founded within the record) exactly Uniform(0, 1).
  const double u_sig  = prior["u_sig"],  u_df   = prior["u_df"];
  const double th_rate = prior["th_rate"];
  const double th_lo = prior["lth_lo"], th_hi = prior["lth_hi"];
  const bool geom = (tstar_p > 0.0);
  const double steps_per_year = (double)M;
  const double z_lo = std::log(ne_lo), z_hi = std::log(ne_hi);

  // log prior on log(nu): half-normal with the Jacobian
  auto lp_lognu = [&](double l) {
    const double nu = std::exp(l);
    return -0.5 * (nu / nu_s) * (nu / nu_s) + l;
  };

  std::vector<double> z(z0.begin(), z0.end());
  double u = u0, lnu = std::log(nu0), mu = mu0, s2 = 0.5;
  // One growth advantage per deme, or one for the panel.  Both use the
  // transform alpha = exp(u/M) - 1; they differ in what the u's are drawn
  // from, a single N(0, tau^2) against a learned N(mu_u, sigma_u^2).
  std::vector<double> uv((std::size_t)D, u0), av((std::size_t)D, 0.0),
                      up((std::size_t)D, 0.0);
  // u_mode 0 = one advantage for the panel
  //        1 = one per deme, drawn from N(mu_u, sigma_u^2) with both
  //            hyperparameters learned -- the same shape as the hierarchy
  //            over the effective sizes, so the demes shrink toward each
  //            other instead of each chasing its own handful of surveys
  const bool u_hier = (u_mode == 1);
  // sigma_u carries a half-t(u_df, u_sig), written as the Huang-Wand scale
  // mixture of two inverse-gammas so both conditionals stay conjugate.  Its
  // mode is at zero, which is the point: the hierarchy has to be able to
  // collapse to a single shared advantage, or comparing it against one asks
  // a question the prior has already answered.
  double mu_u = 0.0, s2_u = 0.25, a_u = u_sig * u_sig;
  // carried on the log scale, like nu
  double lth = geom ? std::log(-std::log1p(-tstar_p) * steps_per_year) : 0.0;
  double p_t = tstar_p;
  auto lp_logth = [&](double l) {         // Exponential prior + Jacobian
    return -th_rate * std::exp(l) + l;
  };
  auto p_of = [&](double l) {             // per-year hazard -> per-step
    return -std::expm1(-std::exp(l) / steps_per_year);
  };
  std::vector<double> L((std::size_t)D), Lp((std::size_t)D),
                      zp((std::size_t)D), buf((std::size_t)T + 1);

  double alpha = std::exp(u / (double)M) - 1.0;
  for (int d = 0; d < D; ++d) av[d] = alpha;
  double nu = std::exp(lnu);
  for (int d = 0; d < D; ++d)
    L[d] = loglik_one(obs[d], nlat_of(z[d], M), av[d], nu * std::exp(z[d]),
                      T, buf, founding, p_t);

  Rcpp::NumericMatrix out(n_iter, 2 * D + 7);
  // `win_*` drives adaptation and is reset at every adaptation step.
  // `acc_*` is the rate reported back, and counts post-warm-up sweeps only:
  // tying it to the adaptation resets made it depend on whether `warm` was a
  // multiple of the adaptation window, and it could then exceed one.
  double acc_z = 0, acc_u = 0, acc_n = 0, acc_t = 0;
  double win_z = 0, win_u = 0, win_n = 0, win_t = 0;
  double step_th = 0.5;
  const int W = 50;
  long guard = 0;

  Rcpp::RNGScope scope;

  for (int it = 0; it < n_iter; ++it) {

    // ---- block 1: the effective sizes ------------------------------------
    double took = 0.0;
    if (shared) {
      // One size for the whole panel: a scalar random walk accepted against
      // the summed likelihood, with the prior that the hierarchy would
      // otherwise put on mu.  There is no between-deme variance to estimate,
      // so block 4 is skipped entirely.
      const double zp0 = z[0] + step_z * R::norm_rand();
      if (zp0 < z_lo || zp0 > z_hi) {
        ++guard;
      } else {
        double sp = 0.0, sc = 0.0;
        const int Np = nlat_of(zp0, M);
        for (int d = 0; d < D; ++d) {
          Lp[d] = loglik_one(obs[d], Np, av[d], nu * std::exp(zp0), T, buf,
                             founding, p_t);
          sp += Lp[d]; sc += L[d];
        }
        const double pri = -0.5 * ((zp0 - mu0) * (zp0 - mu0)
                                 - (z[0] - mu0) * (z[0] - mu0)) / (s0 * s0);
        if (std::log(R::unif_rand()) < sp - sc + pri) {
          for (int d = 0; d < D; ++d) { z[d] = zp0; L[d] = Lp[d]; }
          took = (double)D;
        }
      }
    } else {
    for (int d = 0; d < D; ++d) {
      zp[d] = z[d] + step_z * R::norm_rand();
      if (zp[d] < z_lo || zp[d] > z_hi) { ++guard; continue; }
      Lp[d] = loglik_one(obs[d], nlat_of(zp[d], M), av[d],
                         nu * std::exp(zp[d]), T, buf, founding, p_t);
      const double hier = -0.5 * ((zp[d] - mu) * (zp[d] - mu)
                                - (z[d]  - mu) * (z[d]  - mu)) / s2;
      if (std::log(R::unif_rand()) < Lp[d] - L[d] + hier) {
        z[d] = zp[d]; L[d] = Lp[d]; took += 1.0;
      }
    }
    }
    win_z += took / (double)D;
    if (it >= warm) acc_z += took / (double)D;

    // ---- block 2: the growth advantage -----------------------------------
    if (u_hier) {
      // conditionally independent given everything else, exactly like the
      // effective sizes, so one vectorised pass does all D
      double tk = 0.0;
      for (int d = 0; d < D; ++d) {
        up[d] = uv[d] + step_u * R::norm_rand();
        if (up[d] < u_lo || up[d] > u_hi) { ++guard; continue; }
        const double ap = std::exp(up[d] / (double)M) - 1.0;
        const double Lq = loglik_one(obs[d], nlat_of(z[d], M), ap,
                                     nu * std::exp(z[d]), T, buf, founding,
                                     p_t);
        // each deme shrinks toward mu_u with spread sigma_u
        if (std::log(R::unif_rand()) <
            Lq - L[d] - 0.5 * ((up[d] - mu_u) * (up[d] - mu_u)
                             - (uv[d] - mu_u) * (uv[d] - mu_u)) / s2_u) {
          uv[d] = up[d]; av[d] = ap; L[d] = Lq; tk += 1.0;
        }
      }
      win_u += tk / (double)D;
      if (it >= warm) acc_u += tk / (double)D;
      double su = 0.0;
      for (int d = 0; d < D; ++d) su += uv[d];
      u = su / (double)D;                 // reported as the panel average
      // mu_u gets the same N(0, tau^2) the shared advantage gets, so as
      // sigma_u -> 0 this model becomes exactly the shared one
      const double prec = 1.0 / (tau * tau) + (double)D / s2_u;
      mu_u = R::rnorm((su / s2_u) / prec, std::sqrt(1.0 / prec));
      double S = 0.0;
      for (int d = 0; d < D; ++d) S += (uv[d] - mu_u) * (uv[d] - mu_u);
      // half-t by parameter expansion: sigma_u^2 | a ~ IG(df/2, df/a) and
      // a ~ IG(1/2, 1/u_sig^2), both conditionals inverse-gamma
      s2_u = 1.0 / R::rgamma(0.5 * (u_df + (double)D),
                             1.0 / (u_df / a_u + 0.5 * S));
      a_u  = 1.0 / R::rgamma(0.5 * (u_df + 1.0),
                             1.0 / (1.0 / (u_sig * u_sig) + u_df / s2_u));
    } else {
      const double uq = u + step_u * R::norm_rand();
      if (uq >= u_lo && uq <= u_hi) {
        const double ap = std::exp(uq / (double)M) - 1.0;
        double sp = 0.0, sc = 0.0;
        for (int d = 0; d < D; ++d) {
          Lp[d] = loglik_one(obs[d], nlat_of(z[d], M), ap,
                             nu * std::exp(z[d]), T, buf, founding, p_t);
          sp += Lp[d]; sc += L[d];
        }
        if (std::log(R::unif_rand()) <
            sp - sc - 0.5 * (uq * uq - u * u) / (tau * tau)) {
          u = uq; alpha = ap;
          for (int d = 0; d < D; ++d) { L[d] = Lp[d]; av[d] = ap; }
          win_u += 1.0;
          if (it >= warm) acc_u += 1.0;
        }
      } else {
        ++guard;
      }
    }

    // ---- block 3: log nu -------------------------------------------------
    {
      const double lp_ = lnu + step_nu * R::norm_rand();
      if (lp_ > lnu_lo && lp_ < lnu_hi) {
        const double nup = std::exp(lp_);
        double sp = 0.0, sc = 0.0;
        for (int d = 0; d < D; ++d) {
          Lp[d] = loglik_one(obs[d], nlat_of(z[d], M), av[d],
                             nup * std::exp(z[d]), T, buf, founding, p_t);
          sp += Lp[d]; sc += L[d];
        }
        if (std::log(R::unif_rand()) <
            sp - sc + lp_lognu(lp_) - lp_lognu(lnu)) {
          lnu = lp_; nu = nup;
          for (int d = 0; d < D; ++d) L[d] = Lp[d];
          win_n += 1.0;
          if (it >= warm) acc_n += 1.0;
        }
      } else {
        ++guard;
      }
    }

    // ---- block 3b: the founding hazard, when t* is geometric -------------
    if (geom) {
      const double lp_ = lth + step_th * R::norm_rand();
      if (lp_ > th_lo && lp_ < th_hi) {
        const double pp = p_of(lp_);
        double sp = 0.0, sc = 0.0;
        for (int d = 0; d < D; ++d) {
          Lp[d] = loglik_one(obs[d], nlat_of(z[d], M), av[d],
                             nu * std::exp(z[d]), T, buf, founding, pp);
          sp += Lp[d]; sc += L[d];
        }
        if (std::log(R::unif_rand()) <
            sp - sc + lp_logth(lp_) - lp_logth(lth)) {
          lth = lp_; p_t = pp;
          for (int d = 0; d < D; ++d) L[d] = Lp[d];
          win_t += 1.0;
          if (it >= warm) acc_t += 1.0;
        }
      } else {
        ++guard;
      }
    }

    // ---- block 4: mu and sigma^2, exact conditionals ---------------------
    if (!shared) {
      double zs = 0.0;
      for (int d = 0; d < D; ++d) zs += z[d];
      const double prec = 1.0 / (s0 * s0) + (double)D / s2;
      mu = R::rnorm((mu0 / (s0 * s0) + zs / s2) / prec, std::sqrt(1.0 / prec));
      double S = 0.0;
      for (int d = 0; d < D; ++d) S += (z[d] - mu) * (z[d] - mu);
      s2 = 1.0 / R::rgamma(ig_a + 0.5 * (double)D,
                           1.0 / (ig_b + 0.5 * S));
    }

    // ---- adaptation, warm-up only ---------------------------------------
    if (it < warm && (it + 1) % W == 0) {
      auto adj = [&](double w) {
        double e = (w / (double)W - target) * 2.0;
        if (e >  0.7) e =  0.7;
        if (e < -0.7) e = -0.7;
        return std::exp(e);
      };
      step_z  *= adj(win_z);
      step_u  *= adj(win_u);
      step_nu *= adj(win_n);
      if (geom) step_th *= adj(win_t);
      win_z = win_u = win_n = win_t = 0.0;
    }

    for (int d = 0; d < D; ++d) out(it, d) = z[d];
    out(it, D)     = u;
    out(it, D + 1) = shared ? z[0] : mu;
    out(it, D + 2) = shared ? NA_REAL : std::log(s2);
    out(it, D + 3) = nu;
    out(it, D + 4) = geom ? std::exp(lth) : NA_REAL;   // theta, per year
    for (int d = 0; d < D; ++d) out(it, D + 5 + d) = uv[d];
    out(it, 2 * D + 5) = u_hier ? mu_u : NA_REAL;
    out(it, 2 * D + 6) = u_hier ? std::sqrt(s2_u) : NA_REAL;

    if ((it & 255) == 0) Rcpp::checkUserInterrupt();
  }

  const double post = (double)(n_iter - warm);   // > 0: warm < n_iter
  return Rcpp::List::create(
    Rcpp::_["draws"]  = out,
    Rcpp::_["accept"] = Rcpp::NumericVector::create(
      Rcpp::_["z"] = acc_z / post, Rcpp::_["u"] = acc_u / post,
      Rcpp::_["nu"] = acc_n / post, Rcpp::_["theta"] = acc_t / post),
    Rcpp::_["steps"]  = Rcpp::NumericVector::create(
      Rcpp::_["z"] = step_z, Rcpp::_["u"] = step_u,
      Rcpp::_["nu"] = step_nu, Rcpp::_["theta"] = step_th),
    Rcpp::_["guard_rejects"] = (double)guard);
}
