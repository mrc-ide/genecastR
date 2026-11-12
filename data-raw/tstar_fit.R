## The two fits the "Introduction time" vignette shows.
##
## Regenerate with:
##
##   Rscript data-raw/tstar_fit.R
##
## As with uganda_fit.R this writes into the SOURCE tree, so reinstall
## afterwards or gc_example_tstar_fit() keeps returning the old object.
## R6 objects serialise their methods, so any change to a GCFit method also
## needs this re-run before the shipped fit shows the new behaviour.
library(GeneCastR)

d <- gc_example_data(M = 16)

## ---- 1. one prior for the panel --------------------------------------
panel_ctl <- gc_control(tstar_min = 2012, tstar_mean = 2017,
                        tstar_learn = "shape", tstar_shape = 2)
pf <- GCFit$new(d, chains = 4, iter = 4000, warmup = 2000, seed = 101,
                cores = 4, lead_in = 6 * 16, control = panel_ctl)
set.seed(17)
pp <- pf$project(horizon = 3, n_draws = 400, n_paths = 400,
                 thresholds = c(0.01, 0.05, 0.10))

## ---- 2. a prior that differs by deme ---------------------------------
## ILLUSTRATIVE ONLY.  These are not surveillance dates or published
## estimates; they stand in for the kind of per-deme knowledge a
## programme might hold -- here, a variant believed to have entered from
## the north, reaching the north-east and the south-west later.
north <- c("Koboko", "Arua", "Lamwo", "Kaabong", "Agago")
mid   <- c("Kole", "Amolatar", "Katakwi", "Mubende", "Tororo", "Jinja")
south <- c("Kanungu", "Rukiga", "Kasese", "Hoima", "Kapchorwa")
tmin <- c(setNames(rep(2012, length(north)), north),
          setNames(rep(2013, length(mid)),   mid),
          setNames(rep(2014, length(south)), south))
tmean <- c(setNames(rep(2016, length(north)), north),
           setNames(rep(2017, length(mid)),   mid),
           setNames(rep(2019, length(south)), south))
deme_ctl <- gc_control(tstar_min = tmin, tstar_mean = tmean,
                       tstar_learn = "shape", tstar_shape = 2)
df <- GCFit$new(d, chains = 4, iter = 4000, warmup = 2000, seed = 101,
                cores = 4, lead_in = 6 * 16, control = deme_ctl)
set.seed(17)
dp <- df$project(horizon = 3, n_draws = 400, n_paths = 400,
                 thresholds = c(0.01, 0.05, 0.10))

print(pf); print(df)
saveRDS(list(panel = list(fit = pf, pr = pp),
             deme  = list(fit = df, pr = dp),
             tmin = tmin, tmean = tmean),
        "inst/extdata/tstar_fit.rds", compress = "xz")
