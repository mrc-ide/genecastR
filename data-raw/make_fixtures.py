"""Reference values from the Python implementation, for the R parity tests.

Run from the repository root with the project's interpreter:

    .venv/bin/python -I GeneCastR/data-raw/make_fixtures.py

Everything written here is deterministic, so the R side can assert equality
to ~1e-10 rather than comparing distributions.

N_lat is written out as an integer rather than left to be recomputed from z
on the R side.  The likelihood is piecewise constant in z, so a parameter
sitting on a half-integer lattice boundary would otherwise make the test
depend on the last bit of a decimal round-trip rather than on the kernels.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "intro_model"))
import numpy as np
from core import transition
from uganda_lik import tables, _emis
from reintro_p import transition_p, lam_of_nu
from scipy.special import logsumexp
import uganda_hier as H
H.C_MODE = "one"
from uganda_hier import nlat, load, prep_all, YEARS

OUT = os.path.join(os.path.dirname(__file__), "..", "tests", "testthat",
                   "fixtures")
os.makedirs(OUT, exist_ok=True)
# 17 significant digits is the usual round-trip width for an IEEE double, but
# R's parser mis-rounds some 17-digit decimals by a couple of ulp on this
# platform.  Twenty round-trips exactly, and the extra width costs nothing.
G = "%.21g"


def w(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(x if isinstance(x, str) else G % x
                             for x in r) + "\n")
    print(f"  {name:22s} {len(rows):7d} rows")


# ---- kernels, small enough to store entrywise -------------------------
rows = []
for N, alpha in ((5, 0.0), (5, 0.25), (12, -0.5), (12, 0.0217), (40, 0.1)):
    P = transition(N, alpha)
    for i in range(N + 1):
        for j in range(N + 1):
            rows.append((str(N), G % alpha, str(i), str(j), P[i, j]))
w("kernel.csv", ["N", "alpha", "i", "j", "value"], rows)

# ---- Poisson rows ------------------------------------------------------
rows = []
for N, lam in ((5, 0.0), (5, 0.3), (12, 2.0), (12, 25.0), (40, 0.2718)):
    r = transition_p(N, 0.0, lam)[0]
    for k in range(N + 1):
        rows.append((str(N), G % lam, str(k), r[k]))
w("poisson_row.csv", ["N", "lambda", "k", "value"], rows)

# ---- emissions ---------------------------------------------------------
rows = []
for N, x, n in ((5, 0, 10), (5, 3, 10), (12, 0, 1), (12, 12, 12),
                (40, 7, 44)):
    e = _emis(N, {1: (x, n)})[1]
    for k in range(N + 1):
        rows.append((str(N), str(x), str(n), str(k), e[k]))
w("emission.csv", ["N", "x", "n", "k", "value"], rows)

# ---- the Ugandan panel, several parameter settings ---------------------
M = 16
T = YEARS * M
demes = load(os.path.join(os.path.dirname(__file__), "..", "..",
                          "intro_model", "data", "Ugandan_data.txt"))
NAMES, OBS, CUT = prep_all(demes, M)
D = len(NAMES)

fit = np.load(os.path.join(os.path.dirname(__file__), "..", "..",
                           "intro_model", "cache", "rpt_uganda_c1.npz"),
              allow_pickle=True)
flat = fit["C"][:, int(fit["warm"]):, :].reshape(-1, D + 4)
med = np.median(flat, 0)

SET = {
    # name            u        nu       z
    "median":    (float(med[D]), float(med[D + 3]), med[:D].copy()),
    "neutral":   (0.0,  0.02,  np.full(D, np.log(15.0))),
    "small":     (-1.5, 0.001, np.full(D, np.log(2.5))),
    "large":     (2.0,  0.2,   np.full(D, np.log(100.0))),
}

srows, lrows, trows = [], [], []
for name, (u, nu, z) in SET.items():
    Nl = nlat(z, M)
    al = float(np.exp(u / M) - 1.0)
    lam = lam_of_nu(nu, np.exp(z))
    for d in range(D):
        srows.append((name, G % u, G % nu, G % al, str(M), str(T),
                      NAMES[d], str(int(Nl[d])), lam[d], z[d]))
    for d in range(D):
        P = transition_p(int(Nl[d]), al, float(lam[d]))
        L = tables([OBS[d]], [CUT[d]], int(Nl[d]), al, T, P=P)[:, 1, 0]
        lrows.append((name, NAMES[d], float(logsumexp(L) - np.log(T))))
        post = np.exp(L - logsumexp(L))
        for t in range(T):
            trows.append((name, NAMES[d], str(t + 1), post[t]))

# z is written out too: N_lat is the ROUNDED lattice, so reconstructing z
# from it changes N_eff by up to half a lattice point, and lambda = nu * N_eff
# with it.  The two are not interchangeable.
w("settings.csv",
  ["setting", "u", "nu", "alpha", "M", "T", "deme", "Nlat", "lambda", "z"],
  srows)
w("loglik.csv", ["setting", "deme", "loglik"], lrows)
w("tstar.csv", ["setting", "deme", "t", "prob"], trows)

# ---- the integer N_lat rounding rule ----------------------------------
rows = []
for M_ in (16, 32):
    for z in (np.log(2.0 / M_), np.log(0.5), 0.0, np.log(1.5),
              np.log(2.5), np.log(15.0), np.log(15.5), np.log(100.0),
              # exact half-integers, where the rounding rule shows
              np.log(0.5 / M_ * 3), np.log(2.5 / M_), np.log(4.5 / M_)):
        rows.append((str(M_), G % z, str(int(nlat(np.array([z]), M_)[0]))))
w("nlat.csv", ["M", "z", "Nlat"], rows)
print("done")

# ---- smoother, P(present) and the projection, one parameter setting ----
from uganda_lik import smooth_one

M = 16
T = YEARS * M
HZ = 3 * M
PRIOR_T = np.full(T, 1.0 / T)
u, nu, z = SET["median"]
Nl = nlat(z, M)
al = float(np.exp(u / M) - 1.0)
lam = lam_of_nu(nu, np.exp(z))
THR = [0.01, 0.05, 0.10]

prows, erows = [], []
for d in range(D):
    nl = int(Nl[d])
    P = transition_p(nl, al, float(lam[d]))
    g, gp = smooth_one(OBS[d], nl, al, T, PRIOR_T, 1, P=P)
    pres = 1.0 - gp - g[:, 0]
    k = np.arange(nl + 1) / nl
    fm = g[:, 1:] @ k[1:]
    for t in range(T):
        prows.append((NAMES[d], str(t + 1), pres[t], fm[t]))
        for j, thr in enumerate(THR):
            erows.append((NAMES[d], str(t + 1), G % thr,
                          float(g[t, 1:][k[1:] > thr].sum())))
    p = g[T - 1].copy()
    for h in range(HZ):
        p = p @ P
        prows.append((NAMES[d], str(T + h + 1), float(p[1:].sum()),
                      float(p[1:] @ k[1:])))
        for j, thr in enumerate(THR):
            erows.append((NAMES[d], str(T + h + 1), G % thr,
                          float(p[1:][k[1:] > thr].sum())))

w("smooth.csv", ["deme", "step", "pres", "fmean"], prows)
w("exceed.csv", ["deme", "step", "threshold", "prob"], erows)
