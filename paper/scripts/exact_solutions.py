"""Exact stationary solutions for the validation of the population and lineage layers.

Beentjes, Perez-Carrasco & Grima (Phys. Rev. E 2020, 101:032403) and Jia & Grima
(iScience 2023, 26:105746) solve this class of model exactly and give different
closed forms for a single lineage and for a population snapshot. This module
computes those solutions for the analytically solvable members simulated by
`paper/scripts/fig9_validation.jl`, and compares them with the simulated
distributions.

Two families are covered.

1. A memoryless (exponentially distributed) interdivision time. The simulator
   advances each cell by `dt`, running the chemical master equation exactly over
   the step, and then divides it. Because the interdivision time is memoryless
   and cells are born on the update grid, the number of steps to division is
   geometric with p = 1 - exp(-dt / T), so the recorded distribution is the
   stationary state of one step,

       lineage      pi = [(1-p) I + p B] exp(dt A) pi
       population   rho pi = [(1-p) I + 2 p B] exp(dt A) pi ,   rho = 1 + p

   where A is the reaction generator and B is binomial thinning at division.
   The population operator carries the factor two because both daughters are
   kept, which is what over-weights recently divided (and so depleted) cells.
   The continuum limit of the first two factorial moments is available in closed
   form and is used as an independent check of the matrix solution.

2. A deterministic cycle with exponential volume growth and gene replication.
   Production is volume-scaled and the only non-linear step is division, so the
   count distribution at every cycle phase is exactly Poisson, with a mean that
   solves a linear recursion over the update steps and is halved at division.
   The lineage and the population then differ only in how cycle phases are
   weighted: uniformly along a lineage, and towards younger cells in a growing
   population. The population weights are obtained by iterating the age map of
   the simulated scheme from the founder ages actually used, so no asymptotic
   argument is needed.

Usage: python paper/scripts/exact_solutions.py
Reads paper/output/fig9_counts.csv, fig9_ages.csv and fig9_settings.csv and
writes paper/output/fig9_validation.csv, fig9_pmf.csv and fig9_agefit.csv.
"""

import math
import os
import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.stats import poisson

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output")
NMAX = 250                     # truncation of the count state space


# --------------------------------------------------------------------------- generators
def generator_constitutive(k, gamma, burst=1, nmax=NMAX):
    """Zero-order production in bursts of `burst` molecules, first-order decay."""
    n = nmax + 1
    A = np.zeros((n, n))
    for m in range(n):
        if m + burst < n:
            A[m + burst, m] += k
            A[m, m] -= k
        if m > 0:
            A[m - 1, m] += gamma * m
            A[m, m] -= gamma * m
    return A


def generator_telegraph(k_on, k_off, k_tx, gamma, nmax=NMAX):
    """Two-state promoter, transcription from the active state, first-order decay.

    States are ordered (off, 0), (off, 1), ... , (on, 0), (on, 1), ...
    """
    n = nmax + 1
    A = np.zeros((2 * n, 2 * n))
    off = lambda m: m
    on = lambda m: n + m
    for m in range(n):
        A[on(m), off(m)] += k_on
        A[off(m), off(m)] -= k_on
        A[off(m), on(m)] += k_off
        A[on(m), on(m)] -= k_off
        if m + 1 < n:
            A[on(m + 1), on(m)] += k_tx
            A[on(m), on(m)] -= k_tx
        if m > 0:
            for s in (off, on):
                A[s(m - 1), s(m)] += gamma * m
                A[s(m), s(m)] -= gamma * m
    return A


def binomial_matrix(nmax=NMAX):
    """B[n, m] = probability that a daughter receives n of the mother's m molecules."""
    n = nmax + 1
    B = np.zeros((n, n))
    logfac = np.concatenate(([0.0], np.cumsum(np.log(np.arange(1, n)))))
    for m in range(n):
        idx = np.arange(m + 1)
        logp = logfac[m] - logfac[idx] - logfac[m - idx] - m * np.log(2.0)
        B[idx, m] = np.exp(logp)
    return B


def leading_eigenvector(M, tol=1e-14, maxiter=200_000):
    """Stationary (or fastest-growing) state of a non-negative matrix, by power iteration."""
    v = np.ones(M.shape[0]) / M.shape[0]
    for _ in range(maxiter):
        w = M @ v
        s = w.sum()
        w /= s
        if np.abs(w - v).sum() < tol:
            return w, s
        v = w
    return v, s


def stationary_step(A, p, population, nmax=NMAX, dt=None):
    """Stationary distribution of one update step: exact kinetics over `dt`, then division."""
    L = expm(dt * A)
    B = binomial_matrix(nmax)
    blocks = A.shape[0] // (nmax + 1)
    if blocks == 2:                                   # promoter state is inherited, not partitioned
        Bfull = np.zeros_like(A)
        Bfull[: nmax + 1, : nmax + 1] = B
        Bfull[nmax + 1:, nmax + 1:] = B
        B = Bfull
    H = (1 - p) * np.eye(A.shape[0]) + (2 if population else 1) * p * B
    v, _ = leading_eigenvector(H @ L)
    if blocks == 2:
        v = v[: nmax + 1] + v[nmax + 1:]
    return v / v.sum()


# --------------------------- closed-form factorial moments (continuum limit, for checking)
def factorial_moments(influx_terms, gamma, lam, population, order=2):
    """c_j of G(u) = sum_j c_j u^j, with u = z - 1, for the memoryless-timer models.

    `influx_terms[i]` is the coefficient of u^i contributed by production, so a burst
    of size b gives binomial(b, i) times the burst rate.
    """
    c = [1.0]
    for j in range(1, order + 1):
        den = gamma * j + (lam * (2 - 2.0 ** (1 - j)) if population else lam * (1 - 2.0 ** -j))
        num = sum(influx_terms[i] * c[j - i] for i in range(1, min(len(influx_terms), j + 1)))
        c.append(num / den)
    return c


def telegraph_moments(k_on, k_off, k_tx, gamma, lam, population, order=2):
    """(a_j, b_j) of the off and on generating functions, solved order by order."""
    a = [k_off / (k_on + k_off)]
    b = [k_on / (k_on + k_off)]
    for j in range(1, order + 1):
        d = gamma * j - (lam * (2.0 ** (1 - j) - 2) if population else lam * (2.0 ** -j - 1))
        M = np.array([[-k_on - d, k_off], [k_on, -k_off - d]])
        rhs = np.array([0.0, -k_tx * b[j - 1]])
        aj, bj = np.linalg.solve(M, rhs)
        a.append(aj)
        b.append(bj)
    return [ai + bi for ai, bi in zip(a, b)]


def moments_from_c(c):
    mean = c[1]
    var = 2 * c[2] + c[1] - c[1] ** 2
    return mean, var


# ------------------------------------------------- deterministic cycle with gene replication
def replication_means(k, gamma, lam, dt, J, rep_fraction):
    """Mean count at each recorded cycle phase, for the scheme the simulator runs.

    Volume is held constant within a step and the gene copy number switches at the end
    of the step that crosses `rep_fraction`, so the production rate is a staircase.
    """
    T = J * dt
    e = np.exp(-gamma * dt)
    rate = np.array([k * np.exp(lam * j * dt) * (2.0 if j * dt >= rep_fraction * T else 1.0)
                     for j in range(J)])
    # m_{j+1} = m_j e + (rate_j / gamma)(1 - e);  m_J = A m_0 + B;  m_0 = B / (2 - A)
    A = e ** J
    B = 0.0
    for j in range(J):
        B = B * e + (rate[j] / gamma) * (1 - e)
    m0 = B / (2 - A)
    m = np.empty(J)
    m[0] = m0
    cur = m0
    for j in range(J - 1):
        cur = cur * e + (rate[j] / gamma) * (1 - e)
        m[j + 1] = cur
    return m


def population_age_weights(J, founder_hist, steps):
    """Age distribution of the simulated population after `steps` update steps.

    One step ages every cell and replaces each cell at the last age by two newborns.
    """
    v = np.asarray(founder_hist, dtype=float)
    v = v / v.sum()
    for _ in range(steps):
        w = np.empty_like(v)
        w[1:] = v[:-1]
        w[0] = 2 * v[-1]
        v = w / w.sum()
    return v


def poisson_mixture(weights, means, nmax=NMAX):
    ns = np.arange(nmax + 1)
    return np.array([np.sum(weights * poisson.pmf(n, means)) for n in ns])


# --------------------------------------------------------------------------- comparison
def ks(emp_counts, pmf):
    n = emp_counts.sum()
    ecdf = np.cumsum(emp_counts) / n
    return float(np.max(np.abs(ecdf - np.cumsum(pmf))))


def moments_of(counts):
    ns = np.arange(len(counts))
    n = counts.sum()
    mean = float((ns * counts).sum() / n)
    var = float((counts * (ns - mean) ** 2).sum() / (n - 1))
    return mean, var, int(n)


def histogram(df, nmax=NMAX):
    h = np.zeros(nmax + 1)
    for n, c in zip(df["n"].to_numpy(), df["cells"].to_numpy()):
        h[int(n)] += c
    return h


def main():
    counts = pd.read_csv(os.path.join(OUT, "fig9_counts.csv"))
    cfg = pd.read_csv(os.path.join(OUT, "fig9_settings.csv")).set_index("key")["value"].astype(float)
    k, gamma, burst = cfg["k_prod"], cfg["gamma"], int(cfg["burst"])
    k_on, k_off, k_tx = cfg["k_on"], cfg["k_off"], cfg["k_tx"]
    T_div, dt0 = cfg["T_div"], cfg["dt"]
    J, lam, rep_f = int(cfg["steps_per_cycle"]), cfg["growth_rate"], cfg["replication_fraction"]
    lam_div = 1.0 / T_div

    def exact_pmf(case, mode, dt):
        pop = mode == "population"
        p = 1 - np.exp(-dt / T_div)
        if case in ("constitutive", "constitutive_dt"):
            return stationary_step(generator_constitutive(k, gamma, 1), p, pop, dt=dt)
        if case == "bursty":
            return stationary_step(generator_constitutive(k / burst, gamma, burst), p, pop, dt=dt)
        if case == "telegraph":
            return stationary_step(generator_telegraph(k_on, k_off, k_tx, gamma), p, pop, dt=dt)
        raise KeyError(case)

    def exact_moments(case, mode):
        """Closed-form continuum moments, independent of the matrix solution."""
        pop = mode == "population"
        if case in ("constitutive", "constitutive_dt"):
            return moments_from_c(factorial_moments([0.0, k], gamma, lam_div, pop))
        if case == "bursty":
            terms = [0.0] + [k / burst * float(math.comb(burst, i)) for i in range(1, burst + 1)]
            return moments_from_c(factorial_moments(terms, gamma, lam_div, pop))
        if case == "telegraph":
            return moments_from_c(telegraph_moments(k_on, k_off, k_tx, gamma, lam_div, pop))
        return (np.nan, np.nan)

    rows, pmf_rows = [], []

    # ---- memoryless-timer models -------------------------------------------------
    cache = {}
    for (case, mode, kern, dt), g in counts[counts.case != "replication"].groupby(
            ["case", "mode", "kernel", "dt"], sort=False):
        key = (case, mode, float(dt))
        if key not in cache:
            cache[key] = exact_pmf(case, mode, float(dt))
        pmf = cache[key]
        other = "population" if mode == "lineage" else "lineage"
        okey = (case, other, float(dt))
        if okey not in cache:
            cache[okey] = exact_pmf(case, other, float(dt))
        pmf_other = cache[okey]
        for rep, gr in g.groupby("replicate", sort=False):
            h = histogram(gr)
            mean, var, n = moments_of(h)
            cm, cv = exact_moments(case, mode)
            ex_mean = float((np.arange(len(pmf)) * pmf).sum())
            ex_var = float(((np.arange(len(pmf)) - ex_mean) ** 2 * pmf).sum())
            rows.append(dict(case=case, mode=mode, kernel=kern, dt=float(dt), replicate=int(rep),
                             n_cells=n, mean_sim=mean, mean_exact=ex_mean, mean_continuum=cm,
                             sd_sim=np.sqrt(var), sd_exact=np.sqrt(ex_var), sd_continuum=np.sqrt(cv),
                             ks=ks(h, pmf), ks_other_mode=ks(h, pmf_other),
                             ks_crit99=1.628 / np.sqrt(n)))
        if case == "constitutive" and kern == "DirectSSA" and float(dt) == dt0:
            for n_, (ps, po) in enumerate(zip(pmf, pmf_other)):
                if ps > 1e-9 or po > 1e-9:
                    pmf_rows.append(dict(case=case, mode=mode, n=n_, exact=ps))

    # ---- deterministic cycle with replication ------------------------------------
    rep_counts = counts[counts.case == "replication"]
    if len(rep_counts):
        ages = pd.read_csv(os.path.join(OUT, "fig9_ages.csv"))
        means = replication_means(k, gamma, lam, dt0, J, rep_f)
        founder = np.ones(J)                                   # founders placed uniformly on the grid
        w_lin = np.ones(J) / J
        pmf_lin = poisson_mixture(w_lin, means)
        steps_pop = int(round(cfg["t_end_repl_population"] / dt0))
        w_pop = population_age_weights(J, founder, steps_pop)
        pmf_pop = poisson_mixture(w_pop, means)
        age_rows = []
        for mode, pmf, pmf_other, w in (("lineage", pmf_lin, pmf_pop, w_lin),
                                        ("population", pmf_pop, pmf_lin, w_pop)):
            g = rep_counts[rep_counts["mode"] == mode]
            for rep, gr in g.groupby("replicate", sort=False):
                h = histogram(gr)
                mean, var, n = moments_of(h)
                ex_mean = float((np.arange(len(pmf)) * pmf).sum())
                ex_var = float(((np.arange(len(pmf)) - ex_mean) ** 2 * pmf).sum())
                rows.append(dict(case="replication", mode=mode, kernel="DirectSSA", dt=dt0,
                                 replicate=int(rep), n_cells=n, mean_sim=mean, mean_exact=ex_mean,
                                 mean_continuum=np.nan, sd_sim=np.sqrt(var), sd_exact=np.sqrt(ex_var),
                                 sd_continuum=np.nan, ks=ks(h, pmf), ks_other_mode=ks(h, pmf_other),
                                 ks_crit99=1.628 / np.sqrt(n)))
            a = ages[ages["mode"] == mode]
            if len(a):
                ah = np.zeros(J)
                for step, c in zip(a["age_step"].to_numpy(), a["cells"].to_numpy()):
                    ah[int(step)] += c
                ah = ah / ah.sum()
                for j in range(J):
                    age_rows.append(dict(mode=mode, age_step=j, phase=j / J,
                                         observed=ah[j], expected=w[j]))
        pd.DataFrame(age_rows).to_csv(os.path.join(OUT, "fig9_agefit.csv"), index=False)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "fig9_validation.csv"), index=False)
    if pmf_rows:
        pd.DataFrame(pmf_rows).to_csv(os.path.join(OUT, "fig9_pmf.csv"), index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
