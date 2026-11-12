## 0 = a single copy at t*, 1 = an unconditional Poisson draw there,
## 2 = a Poisson draw conditioned on being non-zero.
.founding_code <- function(x)
  switch(x, one = 0L, rate = 1L, conditional = 2L,
         stop("unknown founding law: ", x, call. = FALSE))

## posterior::subset_draws() reads a bare name as the stem of an indexed
## family, so `u` also hands back every `u[deme]` and `N_eff` every
## `N_eff[deme]`.  Every selection in this package is of exact names, so
## anchor them and escape the index brackets.
.exact_draws <- function(draws, vars) {
  posterior::subset_draws(
    draws, variable = paste0("^", gsub("([][()])", "\\\\\\1", vars), "$"),
    regex = TRUE)
}

## What to call the panel-level effective size.  Under a hierarchy it is the
## exponential of the mean *log* size, so it describes a typical deme rather
## than any particular one, and `exp_mu` says exactly that.  With one shared
## size there is nothing to average over -- the column simply carries that
## size -- so `exp_mu` would be a name borrowed from a parameterisation the
## fit does not use.
.neff_var <- function(control)
  if (identical(control$neff, "shared")) "N_eff" else "exp_mu"
.lneff_var <- function(control)
  if (identical(control$neff, "shared")) "log_N_eff" else "mu"

#' Prior and sampler settings
#'
#' @param tau Scale of the `N(0, tau^2)` prior on `u = log(1 + A)`. At
#'   `tau = 1` about 90% of the prior mass on `A` lies between a fifth and
#'   five times. Neutrality is an interior point with positive density, not a
#'   boundary.
#' @param u_bounds Hard bounds on `u`, applied as a rejection.
#' @param nu_scale Scale of the `HalfNormal(nu_scale)` prior on the influx.
#'   No lower bound is needed: the likelihood goes to `-Inf` as `nu -> 0`,
#'   because with no imports an empty deme cannot produce a non-zero count.
#' @param lognu_bounds Hard bounds on `log(nu)`.
#' @param mu_mean,mu_sd Normal prior on the mean log effective size.
#' @param sigma2_shape,sigma2_rate Inverse-gamma prior on the between-deme
#'   variance. The default `InvGamma(2, 1/2)` puts `P(sigma < 0.2)` at
#'   `5e-5`, which keeps the centred hierarchy out of its funnel. It cannot
#'   express complete sharing, so it is the wrong prior if sharing is the
#'   question being asked.
#' @param u_sigma_scale,u_sigma_df Half-`t` prior on `sigma_u`, the spread of
#'   the per-deme growth advantages under `selection = "hierarchical"`. The
#'   default `u_sigma_df = 1` makes it a half-Cauchy of scale
#'   `u_sigma_scale`, whose **mode is at zero**: the hierarchy can collapse
#'   to a single shared advantage, so comparing the two is a question the
#'   prior has not already answered. An inverse-gamma on `sigma_u^2` cannot
#'   do this -- its density vanishes at zero for any shape and rate -- which
#'   is why this one is not `sigma2_shape`/`sigma2_rate`. Sampled by the
#'   Huang-Wand expansion, so both conditionals stay conjugate and the draw
#'   is exact.
#' @param neff_bounds Hard bounds on `N_eff`, applied **before** the
#'   likelihood is evaluated. They stop a wild proposal allocating a vast
#'   lattice, and the rejections they generate are the negative feedback that
#'   keeps step-size adaptation stable.
#' @param founding What `t*` means. `"conditional"`, the default, makes it
#'   the arrival, with the founding size drawn from `Poisson(lambda)`
#'   conditioned on being non-zero -- the same law that governs every later
#'   re-entry, so the first arrival is not a special case. `"one"` makes it
#'   the arrival with exactly one copy and the flow process
#'   governs only what follows. `"rate"` makes it the step from which the
#'   arrival rate turns positive: the deme joins the flow process at `t*`
#'   and its first actual arrival happens at or after it, drawn from the
#'   same `Poisson(lambda)` that governs every later re-entry. Under
#'   `"one"` and `"conditional"` the first presence time is `t*`
#'   identically; under `"rate"` it is not. `"conditional"` sits between them: `t*` is the arrival, as
#'   under `"one"`, but its size is `Poisson(lambda)` conditioned on being
#'   non-zero rather than forced to a single copy, so the founding event is
#'   drawn from the same law as every later re-entry.
#' @param neff `"shared"` (the default) fits a single effective size for the
#'   whole panel. `sigma` is then not a parameter and is reported as `NA`,
#'   and that single size is reported as `N_eff`, with its log as
#'   `log_N_eff`. `"hierarchical"` gives each deme
#'   its own instead, shrunk toward a common `mu` with spread `sigma`. The
#'   hierarchical prior on `sigma` cannot express complete sharing, so
#'   `"shared"` is the way to ask whether the demes differ at all. It is the
#'   default because deme sizes shrunk hard toward `mu` often buy barely one
#'   effective parameter between them; `gc_waic_compare()` will say whether
#'   yours buy more.
#' @param selection `"shared"` (the default) fits one growth advantage for
#'   the panel. `"hierarchical"` gives each deme its own, drawn from a
#'   learned `N(mu_u, sigma_u^2)`, so the demes shrink toward each other and
#'   can collapse back onto a single advantage. Both share the transform
#'   `alpha = exp(u/M) - 1`. Under the hierarchy the reported `u` and `A` are
#'   panel averages and the per-deme values are `u[deme]` and `A[deme]`.
#' @param tstar Prior on the introduction time. `"geometric"`, the default,
#'   gives founding a constant hazard `theta` per unit time, so `P(t* > T)`
#'   is positive and learned: a deme can simply not have been reached yet.
#'   `"uniform"` spreads `t*` over the observation window instead, under
#'   which the deme is founded inside the record with probability one -- a
#'   deme never seen positive must then have been reached and then hidden or
#'   extinguished, which is rarely what anyone means.
#' @param tstar_rate Rate of the `Exponential` prior on `theta`. The default,
#'   `NULL`, uses the span of the observation window in units of the `time`
#'   column, which makes
#'   `P(founded within the record)` exactly `Uniform(0, 1)` -- no opinion
#'   either way, which is the point of asking. The fit reports `theta` and,
#'   as a transform of it, `p_beyond_baseline = exp(-theta * years)`: the
#'   chance a deme with no data at all has not been reached. Per-deme
#'   posteriors, which use each deme's surveys, are a different quantity and
#'   live in `GCProjection$demes$p_tstar_beyond`.
#' @param lntheta_bounds Hard bounds on `log(theta)`.
#' @param target Acceptance rate the warm-up adaptation aims at.
#' @param step_z,step_u,step_nu Starting random-walk scales.
#' @return A list of settings, for [GCFit].
#' @export
gc_control <- function(tau = 1,
                       founding = c("conditional", "one", "rate"),
                       neff = c("shared", "hierarchical"),
                       selection = c("shared", "hierarchical"),
                       tstar = c("geometric", "uniform"),
                       tstar_rate = NULL, lntheta_bounds = c(-8, 3),
                       u_bounds = c(-4, 4),
                       nu_scale = 0.02,
                       lognu_bounds = c(-14, 2),
                       mu_mean = log(15), mu_sd = 1.5,
                       sigma2_shape = 2, sigma2_rate = 0.5,
                       u_sigma_scale = 0.5, u_sigma_df = 1,
                       neff_bounds = c(2, 400),
                       target = 0.30,
                       step_z = 0.9, step_u = 0.25, step_nu = 0.6) {
  founding <- match.arg(founding)
  tstar <- match.arg(tstar)
  neff <- match.arg(neff)
  selection <- match.arg(selection)
  stopifnot(tau > 0, nu_scale > 0, mu_sd > 0, sigma2_shape > 0,
            u_sigma_scale > 0, u_sigma_df > 0,
            sigma2_rate > 0, target > 0, target < 1,
            length(u_bounds) == 2L, u_bounds[1L] < u_bounds[2L],
            length(neff_bounds) == 2L, neff_bounds[1L] >= 2,
            neff_bounds[1L] < neff_bounds[2L])
  list(tau = tau, u_bounds = u_bounds,
       nu_scale = nu_scale, lognu_bounds = lognu_bounds,
       mu_mean = mu_mean, mu_sd = mu_sd,
       sigma2_shape = sigma2_shape, sigma2_rate = sigma2_rate,
       u_sigma_scale = u_sigma_scale, u_sigma_df = u_sigma_df,
       neff_bounds = neff_bounds,
       founding = founding, neff = neff, selection = selection,
       tstar = tstar,
       tstar_rate = tstar_rate,
       lntheta_bounds = lntheta_bounds,
       target = target,
       step_z = step_z, step_u = step_u, step_nu = step_nu)
}

## `gc_control()` is user-facing, so it names its settings as its arguments do
## and keeps each pair of bounds in one vector.  The compiled sampler wants
## flat scalars under its own shorter names; that translation lives here, in
## one place, rather than leaking into what the user reads back off a fit.
.prior_list <- function(control, th_rate) {
  list(tau    = control$tau,
       u_lo   = control$u_bounds[1L],       u_hi   = control$u_bounds[2L],
       nu_s   = control$nu_scale,
       lnu_lo = control$lognu_bounds[1L],   lnu_hi = control$lognu_bounds[2L],
       mu0    = control$mu_mean,            s0     = control$mu_sd,
       ig_a   = control$sigma2_shape,       ig_b   = control$sigma2_rate,
       ne_lo  = control$neff_bounds[1L],    ne_hi  = control$neff_bounds[2L],
       u_sig  = control$u_sigma_scale,      u_df   = control$u_sigma_df,
       lth_lo = control$lntheta_bounds[1L], lth_hi = control$lntheta_bounds[2L],
       th_rate = th_rate)
}


#' A fitted Wright-Fisher introduction model
#'
#' Runs blocked Metropolis-within-Gibbs against the exact likelihood. Because
#' the likelihood is exact rather than estimated, this is ordinary
#' Metropolis-within-Gibbs and not a pseudo-marginal scheme: the acceptance
#' ratio uses the true value.
#'
#' Chains start dispersed on purpose, over a fortyfold range in the influx
#' and across the sign of `u`, so that the convergence diagnostic has
#' something to detect.
#'
#' @examples
#' \donttest{
#' d <- gc_example_data()
#' fit <- GCFit$new(d, chains = 2, iter = 400)
#' fit$summary()
#' }
#'
#' @export
GCFit <- R6::R6Class(
  "GCFit",
  public = list(

    #' @field data The [GCData] the model was fitted to.
    data = NULL,
    #' @field draws A `posterior::draws_array`: iterations x chains x
    #'   variables, warm-up already discarded.
    draws = NULL,
    #' @field control The settings used, from [gc_control].
    control = NULL,
    #' @field info Per-chain acceptance rates, final step sizes, guard
    #'   rejections and run time.
    info = NULL,

    #' @description Fit the model.
    #' @param data A [GCData].
    #' @param chains Number of chains.
    #' @param iter Sweeps per chain, warm-up included.
    #' @param warmup Warm-up sweeps; defaults to a third of `iter`. Step
    #'   sizes adapt only within these, and the draws are discarded.
    #' @param seed Base seed; chain `i` uses `seed + i`.
    #' @param control Settings from [gc_control].
    #' @param cores Chains to run at once. Chains are independent, so this
    #'   is the cheapest speedup available. The default of 1 is deliberate:
    #'   package checks must not take more than two cores, and a caller who
    #'   wants more should say so.
    #' @param lead_in Lead-in before the first survey, **in steps**. This
    #'   is where `t0` is placed, and so the earliest arrival the model will
    #'   entertain: a three-year lead-in at `M = 16` is `lead_in = 3 * 16`.
    #'   The default, `NULL`, keeps whatever window `data` was built with.
    #'   Any other value re-anchors the panel, rebuilding the [GCData]; the
    #'   window actually fitted is always `fit$data`.
    #' @param quiet Suppress progress messages.
    initialize = function(data, chains = 4L, iter = 3000L, warmup = NULL,
                          seed = 1L, control = gc_control(), cores = 1L,
                          lead_in = NULL, quiet = FALSE) {
      if (!inherits(data, "GCData")) stop("`data` must be a GCData", call. = FALSE)
      if (!is.null(lead_in)) data <- private$reanchor(data, lead_in, quiet)
      chains <- as.integer(chains); iter <- as.integer(iter)
      if (chains < 1L) stop("`chains` must be at least 1", call. = FALSE)
      if (iter < 10L) stop("`iter` must be at least 10", call. = FALSE)
      warmup <- if (is.null(warmup)) iter %/% 3L else as.integer(warmup)
      if (warmup < 0L || warmup >= iter)
        stop("`warmup` must be in [0, iter)", call. = FALSE)

      D <- data$n_demes
      dd <- data$deme_data()
      ## Exponential(years) on theta makes P(founded within the record)
      ## uniform, which is the honest default when the question is whether
      ## it has arrived at all.
      ## the window's span in whatever unit `time` is measured in; theta is
      ## a hazard per that same unit
      span <- data$span
      th_rate <- if (is.null(control$tstar_rate)) span else control$tstar_rate
      prior <- .prior_list(control, th_rate)
      geom <- identical(control$tstar, "geometric")
      ## start at the prior mean hazard, converted to a per-step probability
      tstar_p0 <- if (geom) -expm1(-(1 / th_rate) / data$M) else -1

      ## dispersed starts: u across the sign, nu over a fortyfold range
      u0  <- rep_len(c(-0.5, 0.2, 0.6, 1.2), chains)
      nu0 <- rep_len(c(0.0005, 0.002, 0.008, 0.02), chains)

      ## Seeding inside the worker, not outside it: a forked child inherits
      ## the parent's RNG stream, so a chain's draws must be fixed by its own
      ## seed if the fit is to be reproducible whatever `cores` is set to.
      one <- function(i) {
        set.seed(seed + i)
        gc_mwg_run(
          dd, M = data$M, T = data$T, n_iter = iter, warm = warmup,
          z0 = rep(control$mu_mean, D), u0 = u0[i], nu0 = nu0[i],
          step_z = control$step_z, step_u = control$step_u,
          step_nu = control$step_nu, prior = prior, target = control$target,
          founding = .founding_code(control$founding),
          tstar_p = tstar_p0,
          shared = identical(control$neff, "shared"),
          u_mode = switch(control$selection, shared = 0L, hierarchical = 1L))
      }

      cores <- max(1L, min(as.integer(cores), chains))
      t0 <- proc.time()[["elapsed"]]
      if (cores > 1L && .Platform$OS.type == "unix") {
        if (!quiet) message("running ", chains, " chains on ", cores,
                            " cores ...")
        res <- parallel::mclapply(seq_len(chains), one, mc.cores = cores)
        bad <- vapply(res, inherits, logical(1L), "try-error")
        if (any(bad))
          stop("chain(s) ", paste(which(bad), collapse = ", "), " failed: ",
               conditionMessage(attr(res[[which(bad)[1L]]], "condition")),
               call. = FALSE)
      } else if (cores > 1L) {
        ## no fork on Windows, so ship the work to a socket cluster
        if (!quiet) message("running ", chains, " chains on ", cores,
                            " workers ...")
        cl <- parallel::makePSOCKcluster(cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterEvalQ(cl, requireNamespace("GeneCastR", quietly = TRUE))
        res <- parallel::parLapply(cl, seq_len(chains), one)
      } else {
        res <- vector("list", chains)
        for (i in seq_len(chains)) {
          if (!quiet) message("chain ", i, "/", chains, " ...")
          res[[i]] <- one(i)
        }
      }
      wall <- proc.time()[["elapsed"]] - t0

      keep <- (warmup + 1L):iter
      uh <- identical(control$selection, "hierarchical")
      perd <- uh
      nex <- (if (geom) 9L else 7L) + (if (perd) 2L * D else 0L) +
             (if (uh) 2L else 0L)
      arr <- array(NA_real_, dim = c(length(keep), chains, D + nex))
      for (i in seq_len(chains)) {
        m <- res[[i]]$draws[keep, , drop = FALSE]
        arr[, i, seq_len(D)] <- m[, seq_len(D)]      # z_d
        arr[, i, D + 1L] <- m[, D + 1L]              # u
        arr[, i, D + 2L] <- m[, D + 2L]              # mu
        arr[, i, D + 3L] <- m[, D + 3L]              # log sigma^2
        arr[, i, D + 4L] <- m[, D + 4L]              # nu
        arr[, i, D + 5L] <- expm1(m[, D + 1L])       # A
        arr[, i, D + 6L] <- exp(m[, D + 3L] / 2)     # sigma
        arr[, i, D + 7L] <- exp(m[, D + 2L])         # exp(mu)
        if (perd) {
          k <- (if (geom) 9L else 7L)
          arr[, i, D + k + seq_len(D)] <- m[, D + 5L + seq_len(D)]
          arr[, i, D + k + D + seq_len(D)] <- expm1(m[, D + 5L + seq_len(D)])
          if (uh) {
            arr[, i, D + k + 2L * D + 1L] <- m[, 2L * D + 6L]   # mu_u
            arr[, i, D + k + 2L * D + 2L] <- m[, 2L * D + 7L]   # sigma_u
          }
        }
        if (geom) {
          arr[, i, D + 8L] <- m[, D + 5L]            # theta, per unit time
          ## A deterministic transform of theta, like A of u: the chance a
          ## deme with NO data has not been reached within the record.  It
          ## is the baseline the per-deme posteriors in
          ## `GCProjection$demes$p_tstar_beyond` update away from, and it
          ## carries no information theta does not.
          arr[, i, D + 9L] <- exp(-m[, D + 5L] * span)
        }
      }
      ## N_eff on the natural scale, not the sampled log scale: every
      ## reported quantity uses it, and the transformation is where sign
      ## errors live.
      arr[, , seq_len(D)] <- exp(arr[, , seq_len(D), drop = FALSE])
      dimnames(arr) <- list(NULL, NULL,
        c(paste0("N_eff[", data$demes, "]"),
          "u", .lneff_var(control), "log_sigma2", "nu", "A", "sigma",
          .neff_var(control),
          if (geom) c("theta", "p_beyond_baseline"),
          if (perd) c(paste0("u[", data$demes, "]"),
                      paste0("A[", data$demes, "]")),
          if (uh) c("mu_u", "sigma_u")))

      self$data <- data
      self$control <- control
      self$draws <- posterior::as_draws_array(arr)
      self$info <- data.frame(
        chain = seq_len(chains),
        accept_z = vapply(res, function(r) r$accept[["z"]], numeric(1L)),
        accept_u = vapply(res, function(r) r$accept[["u"]], numeric(1L)),
        accept_nu = vapply(res, function(r) r$accept[["nu"]], numeric(1L)),
        accept_theta = vapply(res, function(r) r$accept[["theta"]],
                              numeric(1L)),
        step_z = vapply(res, function(r) r$steps[["z"]], numeric(1L)),
        step_u = vapply(res, function(r) r$steps[["u"]], numeric(1L)),
        step_nu = vapply(res, function(r) r$steps[["nu"]], numeric(1L)),
        step_theta = vapply(res, function(r) r$steps[["theta"]], numeric(1L)),
        guard_rejects = vapply(res, function(r) r$guard_rejects, numeric(1L)),
        row.names = NULL
      )
      attr(self$info, "wall") <- wall
      attr(self$info, "cores") <- cores
      attr(self$info, "iter") <- iter
      attr(self$info, "warmup") <- warmup
      invisible(self)
    },

    #' @description Posterior summary with convergence diagnostics.
    #' @param vars Variable names; defaults to the shared parameters, each
    #'   on the scale it is reported on rather than the scale it is sampled
    #'   on. Pass `NULL` for every variable, `"mu"` and `"log_sigma2"`
    #'   included.
    #' @return A data frame, one row per variable, from the **posterior**
    #'   package -- split-Rhat and both effective sample sizes, not
    #'   hand-rolled.
    summary = function(vars = NULL) {
      if (is.null(vars)) {
        vars <- c("u", "A", "nu", .neff_var(self$control), "sigma")
        if (identical(self$control$selection, "hierarchical"))
          vars <- c(vars, "mu_u", "sigma_u")
        if (identical(self$control$neff, "shared"))
          vars <- setdiff(vars, "sigma")
        if (identical(self$control$tstar, "geometric"))
          vars <- c(vars, "theta", "p_beyond_baseline")
      }
      if (!requireNamespace("posterior", quietly = TRUE))
        stop("the posterior package is needed for summaries", call. = FALSE)
      d <- if (is.null(vars)) self$draws else
        .exact_draws(self$draws, vars)
      as.data.frame(posterior::summarise_draws(d))
    },

    #' @description Smoothed curves, projections and the exceedance table.
    #'
    #' Memoised: a second call with the same arguments returns the stored
    #' object rather than recomputing it. That is why this class has
    #' reference semantics.
    #'
    #' @param horizon Time units to project beyond the last step.
    #' @param n_draws Posterior draws to average over.
    #' @param thresholds Frequencies at which to report exceedance.
    #' @param probs Quantile levels for the frequency columns.
    #' @param nbins Bins for the frequency density.
    #' @param n_paths Paths per deme per draw for `T_est`. Zero (the
    #'   default) skips the path sampler: `T_first` is analytic and costs
    #'   nothing, `T_est` is the only expensive output here.
    #' @return A [GCProjection].
    project = function(horizon = 3, n_draws = 120L,
                       thresholds = c(0.01, 0.05, 0.10),
                       probs = c(.025, .25, .5, .75, .975),
                       nbins = 250L, n_paths = 0L) {
      key <- paste(horizon, n_draws, nbins, n_paths,
                   paste(sort(thresholds), collapse = ","),
                   paste(sort(probs), collapse = ","), sep = "|")
      if (!is.null(private$.proj[[key]])) return(private$.proj[[key]])
      p <- GCProjection$new(self, horizon = horizon, n_draws = n_draws,
                            thresholds = thresholds, probs = probs,
                            nbins = nbins, n_paths = n_paths)
      private$.proj[[key]] <- p
      p
    },

    #' @description Drop everything memoised by `$project()`.
    clear_cache = function() {
      private$.proj <- list()
      invisible(self)
    },

    #' @description Posterior over the introduction time.
    #'
    #' Analytic, not sampled: each draw's contribution is a normalisation of
    #' the backward row the likelihood already sums, so this carries no
    #' Monte-Carlo error beyond the posterior draws themselves.
    #'
    #' @param n_draws Draws to average over, evenly spaced through the chains.
    #' @return A data frame with `deme`, `step`, `time` and `prob`, one row
    #'   per deme and step inside the record. Its `"p_beyond"` attribute
    #'   holds, per deme, the probability that `t*` postdates the record --
    #'   zero under the uniform prior, and zero for any deme ever seen
    #'   positive. `prob` and `p_beyond` together sum to one.
    tstar = function(n_draws = 200L) {
      p <- private$par_draws(n_draws)
      D <- self$data$n_demes; T <- self$data$T
      acc <- matrix(0, D, T + 1L)
      for (i in seq_len(nrow(p$z))) {
        acc <- acc + gc_tstar_pt(gc_nlat(p$z[i, ], self$data$M),
                                 gc_alpha(p$u[i], self$data$M),
                                 gc_lambda(p$nu[i], exp(p$z[i, ])),
                                 T, private$dd(),
                                 .founding_code(self$control$founding),
                                 p$tp[i])
      }
      acc <- acc / nrow(p$z)
      beyond <- stats::setNames(acc[, T + 1L], self$data$demes)
      out <- data.frame(
        deme = rep(self$data$demes, times = T),
        step = rep(seq_len(T), each = D),
        time = rep(self$data$t0 + seq_len(T) / self$data$M, each = D),
        prob = as.vector(acc[, seq_len(T), drop = FALSE]),
        stringsAsFactors = FALSE
      )
      ## the rows cover t* inside the record; what is left over is the
      ## probability the variant had not reached the deme at all, which is
      ## zero under the uniform prior
      attr(out, "p_beyond") <- beyond
      out
    },

    #' @description Print a compact description.
    #' @param ... Ignored.
    print = function(...) {
      nd <- posterior::ndraws(self$draws)
      cat("<GCFit>  Poisson flow, t* ",
          if (identical(self$control$tstar, "geometric"))
            "geometric (may postdate the record)" else "uniform on the window",
          "\n          (",
          switch(self$control$founding,
                 one = "t* = the arrival, one copy",
                 rate = "t* = onset of the arrival rate",
                 conditional = "t* = the first non-zero draw"),
          ")\n", sep = "")
      cat("  data    : ", self$data$n_demes, " demes, ",
          self$data$n_surveys, " surveys, M = ", self$data$M,
          ", T = ", self$data$T, "\n", sep = "")
      cat("  window  : ", format(self$data$t0), " to ",
          format(self$data$t0 + self$data$span), ", lead-in ",
          round(self$data$M * (min(self$data$surveys$time) - self$data$t0)),
          " steps\n", sep = "")
      cat("  chains  : ", posterior::nchains(self$draws), " x ",
          attr(self$info, "iter"), " (", attr(self$info, "warmup"),
          " warm-up), ", nd, " draws retained\n", sep = "")
      cat("  time    : ", sprintf("%.1f s", attr(self$info, "wall")),
          if (attr(self$info, "cores") > 1L)
            paste0(" on ", attr(self$info, "cores"), " cores") else "",
          "\n", sep = "")
      s <- self$summary()
      cat("\n")
      print(s[, c("variable", "median", "q5", "q95", "rhat", "ess_bulk")],
            row.names = FALSE, digits = 3)
      bad <- s$rhat > 1.05
      if (any(bad, na.rm = TRUE)) {
        cat("\n  Rhat above 1.05 for: ",
            paste(s$variable[bad], collapse = ", "),
            " -- run longer.\n", sep = "")
      }
      invisible(self)
    }
  ),

  private = list(
    .proj = list(),
    ## `lead_in` counts steps, not time units, and a whole number of them
    ## keeps every survey on a step whatever window the panel arrived with:
    ## each step index is just the old one shifted by a constant integer.
    reanchor = function(data, lead_in, quiet) {
      if (!(length(lead_in) == 1L && is.finite(lead_in) &&
            lead_in == round(lead_in) && lead_in >= 1)) {
        stop("`lead_in` is a whole number of steps >= 1; a three-year ",
             "lead-in at M = ", data$M, " is lead_in = ", 3L * data$M, ".",
             call. = FALSE)
      }
      lead_in <- as.integer(lead_in)
      t_first <- min(data$surveys$time)
      origin <- t_first - lead_in / data$M
      if (!quiet && !isTRUE(all.equal(origin, data$t0))) {
        message("lead-in of ", lead_in, " steps: t0 moved from ",
                format(data$t0), " to ", format(origin), ".")
      }
      GCData$new(data$surveys, M = data$M, t0 = origin, quiet = TRUE)
    },
    dd = function() {
      if (is.null(private$.dd)) private$.dd <- self$data$deme_data()
      private$.dd
    },
    .dd = NULL,
    par_draws = function(n_draws) {
      m <- posterior::as_draws_matrix(self$draws)
      idx <- unique(round(seq(1, nrow(m), length.out = n_draws)))
      D <- self$data$n_demes
      geom <- identical(self$control$tstar, "geometric")
      list(z  = log(m[idx, seq_len(D), drop = FALSE]),
           u  = as.numeric(m[idx, "u"]),
           nu = as.numeric(m[idx, "nu"]),
           ## each draw's own hazard: collapsing these to the median first
           ## is not the same average, and disagreed with GCProjection
           tp = if (geom) -expm1(-as.numeric(m[idx, "theta"]) / self$data$M)
                else rep(-1, length(idx)))
    }
  )
)
