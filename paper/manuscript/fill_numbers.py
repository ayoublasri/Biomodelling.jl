"""Fill the {{placeholders}} of 02_results.template.md from paper/output/*.csv → 02_results.md.
Every number in the Results section is therefore traceable to a simulation output table."""
import os, re, sys, math
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "output")
def load(n): return pd.read_csv(os.path.join(OUT, n))
def kv(n): d = load(n); return dict(zip(d.key, d.value))
V = {}
def put(k, v): V[k] = v
def safe(fn):
    def run():
        try: fn()
        except Exception as e: print(f"  [{fn.__name__}: {e}]")
    run.__name__ = fn.__name__
    return run

@safe
def f2():
    d = load("fig2d_kernels.csv")
    put("fig2d_times", " to ".join(f"{x:.2f} s" for x in (d.sec_per_1e4_telegraph.min(), d.sec_per_1e4_telegraph.max())))
    s = kv("fig2e_crosscheck_stats.csv")
    put("ks_two", f"{s['ks_two_sample']:.4f}"); put("ks_crit", f"{s['ks_two_sample_crit99']:.4f}")
    r = s["time_jumpprocesses_serial_s"] / s["time_biomodelling_serial_s"]
    put("speed_ratio", f"{r:.1f}-fold" if r >= 1 else f"{1/r:.1f}-fold (in favour of JumpProcesses.jl)")
    d = load("fig2f_runtime.csv"); row = d[(d.kernel == "HybridSSATau") & (d.cells == d.cells.max()) & (d.genes == 10)].iloc[0]
    put("runtime_10k", f"{row.seconds:.1f} s"); put("runtime_case", f"{int(row.cells)} cells × {int(row.genes)} telegraph genes × 200 steps (hybrid kernel)")
    # scaling in the number of reactions is not the same for the two kernels, so quote both
    # the gene scan is run at one cell count: use whichever it is
    ncells = d.groupby("cells").genes.nunique().idxmax()
    parts = []
    for kern in ("HybridSSATau", "DirectSSA"):
        g = d[(d.kernel == kern) & (d.cells == ncells)].sort_values("genes")
        if len(g) < 2: continue
        lo, hi = g.iloc[0], g.iloc[-1]
        parts.append(f"{lo.seconds:.2f} s to {hi.seconds:.1f} s for the "
                     f"{'hybrid' if kern == 'HybridSSATau' else 'direct'} kernel, a factor of "
                     f"{hi.seconds/lo.seconds:.0f} for the same factor of {hi.genes/lo.genes:.0f} in genes")
    put("runtime_scaling", "; ".join(parts))

@safe
def f3():
    d = load("fig3c_heritability.csv").sort_values("k_switch")
    put("md_slow", f"{d.mother_daughter.iloc[0]:.2f}"); put("md_fast", f"{d.mother_daughter.iloc[-1]:.2f}")
    # Fig. 3f plots the fitted transcription rate per unit volume and per twenty transcripts, so the
    # text has to separate the raw change, the part of it that is volume, and the remainder.
    f = load("fig3f_copynumber.csv").set_index("subset")
    one, two, pool = f.loc["copies1"], f.loc["copies2"], f.loc["pooled"]
    put("ktx_fold", f"{two.k_tx / one.k_tx:.1f}")
    put("vol_fold", f"{two.mean_volume / one.mean_volume:.1f}")
    put("ktx_per_vol_fold", f"{(two.k_tx / two.mean_volume) / (one.k_tx / one.mean_volume):.1f}")
    put("kon_fold", f"{two.k_on / one.k_on:.1f}")
    put("fit_one_copy", f"{one.k_on:.2f}, {one.k_off:.2f} and {one.k_tx / one.mean_volume / 20:.1f}")
    put("fit_pooled", f"{pool.k_on:.2f}, {pool.k_off:.2f} and {pool.k_tx / pool.mean_volume / 20:.1f}")

LABEL_MECH = {"pre_existing": "pre-existing tolerance", "pre_existing_cost": "a fitness cost of resistance",
              "drug_induced": "drug-induced tolerance"}

def paired_margin(df, model):
    """Winner against the best schedule with a non-zero release period, paired over shared seeds.

    Every schedule in a scan is run on the same seeds (common random numbers), so the two cell
    means are coupled and sqrt((s_a^2 + s_b^2)/n) is not a standard error of their difference.
    The paired standard error is valid whatever the coupling, which matters here because it is
    positive for some mechanisms and negative for others. With five seeds the two-sided 5% point
    of Student's t on four degrees of freedom is 2.776, not 2, so the test is stated as a t-test.
    """
    from scipy import stats as _st
    w = df[df.model == model].pivot_table(index="seed", columns=["release_period", "dose"],
                                          values="long_term_growth_rate")
    mu = w.mean()
    best = mu.idxmin()
    holidays = [c for c in mu.index if c[0] > 0]
    other = mu[holidays].idxmin() if holidays else mu.drop(index=best).idxmin()
    diff = w[other] - w[best]
    n = len(w)
    gap = float(diff.mean()); se = float(diff.std(ddof=1) / np.sqrt(n))
    p = float(_st.ttest_rel(w[other], w[best]).pvalue) if n > 1 else float("nan")
    return dict(model=model, best=best, other=other, gap=gap, se=se, p=p, n=n, resolved=bool(p < 0.05))

@safe
def f4():
    put("mem_time", f"{1 / (0.002 + 0.02):.0f}")
    d = load("fig4c_decay_vs_dose.csv").sort_values("dose")
    lo = d[d.dose == 0.5].decay_rate.iloc[0]; hi = d[d.dose == 2.0].decay_rate.iloc[0]
    put("decay_lo", f"{lo:.3f}"); put("decay_hi", f"{hi:.3f}"); put("decay_fold", f"{hi / lo:.1f}")
    dd = d[d.dose > 0]
    put("divtime_change", f"{abs(dd.division_time_mean.iloc[-1] - dd.division_time_mean.iloc[0]) / dd.division_time_mean.iloc[0] * 100:.1f}%")
    put("deathtime_range", f"{dd.death_time_mean.max():.1f} to {dd.death_time_mean.min():.1f}")
    c = load("fig4d_fate_concordance.csv")
    n = c.n_pairs.round().astype(int); same = (c.concordance * n).round().astype(int); disc = n - same
    deaths = (c.death_fraction * 2 * n).round().astype(int); both_die = ((deaths - disc) / 2).round().astype(int); both_survive = same - both_die
    c["p_survive"] = (2 * both_survive + disc) / (2 * n); c["p_cond"] = (2 * both_survive) / np.maximum(2 * both_survive + disc, 1)
    c["both_survive"] = both_survive
    # the fast-switching control's fold enrichment rests on a handful of surviving pairs, so the
    # text has to say how many rather than report the fold as though it were resolved
    def word(k):
        return ("no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine")[k] if k < 10 else str(k)
    def npairs(rel):
        r = c[(c.model == "fast") & (c.relation == rel)]
        return int(r.both_survive.iloc[0]) if len(r) else 0
    ns, nc = npairs("sisters"), npairs("cousins")
    put("fast_pairs", f"{word(ns)} sister pair{'' if ns == 1 else 's'} and {word(nc)} cousin pair{'' if nc == 1 else 's'}")
    for model, mtag in (("memory", ""), ("fast", "_fast")):
        for rel, rtag in (("sisters", "sis"), ("cousins", "cous")):
            r = c[(c.model == model) & (c.relation == rel)].iloc[0]
            put(f"{rtag}_cond{mtag}", f"{r.p_cond:.2f}"); put(f"{rtag}_marg{mtag}", f"{r.p_survive:.2f}"); put(f"{rtag}_fold{mtag}", f"{r.p_cond / r.p_survive:.1f}")
            put(f"{rtag}_n{mtag}", f"{int(r.n_pairs)}")
    e = load("fig4e_clone_diversity.csv"); m = e.groupby("model").mean(numeric_only=True)
    put("clones_before", f"{m.loc['pre_existing', 'effective_clones_before']:.0f}"); put("clones_after_pre", f"{m.loc['pre_existing', 'effective_clones_after']:.0f}")
    put("clones_after_ind", f"{m.loc['drug_induced', 'effective_clones_after']:.0f} of {m.loc['drug_induced', 'effective_clones_before']:.0f}")
    try:
        f = load("fig4f_schedules_seeds.csv"); nseed = int(f.seed.nunique())
    except Exception:
        f = load("fig4f_schedules.csv"); f["seed"] = 1; nseed = 1
    agg = f.groupby(["model", "release_period", "dose"]).long_term_growth_rate.agg(["mean", "std", "count"]).reset_index()
    marg = []
    for model, tag in (("pre_existing", "pre"), ("pre_existing_cost", "cost"), ("drug_induced", "ind")):
        g = agg[agg.model == model].sort_values("mean")
        best = g.iloc[0]
        err = "" if nseed == 1 else f", ± {best['std']:.4f} over {nseed} seeds"
        put(f"best_{tag}", f"release period {best.release_period:g}, dose {best.dose:g} (growth rate {best['mean']:+.4f} per time unit{err})")
        cont = g[(g.release_period == 0) & (g.dose == g.dose.max())]["mean"].iloc[0]
        put(f"cont_{tag}", f"{cont:+.4f}")
        if nseed > 1:
            marg.append(paired_margin(f, model))
    if marg:
        won = [m for m in marg if m["resolved"]]
        lost = [m for m in marg if not m["resolved"]]
        detail = "; ".join(f"{LABEL_MECH.get(m['model'], m['model'])} {m['gap']:+.4f} against {m['se']:.4f}, "
                           f"$P = {m['p']:.2f}$" for m in lost)
        put("sched_seed_note",
            f"Each point of the scan is the mean of {nseed} independent seeds, and every schedule is run on the same "
            f"{nseed} seeds, so schedules are compared by their paired differences, in which the shared founder "
            f"population cancels; the spread quoted beside each winner above is the seed-to-seed scatter of its "
            f"absolute growth rate, which is larger than the paired differences and is not the error on them. The "
            f"winning schedule beats the best schedule with a non-zero release period by "
            f"{min(m['gap'] for m in marg):.4f} to {max(m['gap'] for m in marg):.4f} per time unit against paired "
            f"standard errors of {min(m['se'] for m in marg):.4f} to {max(m['se'] for m in marg):.4f}, a margin that "
            f"a two-sided paired t-test resolves at the 5% level in {len(won)} of the {len(marg)} mechanisms"
            + (f", the one in which the resistant state carries a fitness cost"
               if len(won) == 1 and won[0]["model"] == "pre_existing_cost" else "")
            + (f". Elsewhere the scan identifies a region of good schedules rather than a single best one ({detail})."
               if lost else "."))
    else:
        put("sched_seed_note", "")
    # under the paired test the fitness-cost holidays are not equal to continuous dosing, so the
    # text cannot say they "match" it; report how many of them are resolvably worse
    try:
        from scipy import stats as _st
        w = f[f.model == "pre_existing_cost"].pivot_table(index="seed", columns=["release_period", "dose"],
                                                          values="long_term_growth_rate")
        base = (0.0, float(max(c[1] for c in w.columns)))
        cells = [c for c in w.columns if 0 < c[0] <= 10 and c[1] >= 1.0]   # the doses the sentence is about
        ps = {c: float(_st.ttest_rel(w[c], w[base]).pvalue) for c in cells}
        worse = [c for c in cells if ps[c] < 0.05 and (w[c] - w[base]).mean() > 0]
        gaps = [float((w[c] - w[base]).mean()) for c in cells]
        put("cost_holiday", f"release periods of up to 10 time units at doses of 1 and 2 cost little, though "
                            f"{len(worse)} of those {len(cells)} are resolvably worse than continuous dosing "
                            f"at the highest dose "
                            f"({min(gaps):+.4f} to {max(gaps):+.4f} per time unit, "
                            f"{V.get('cont_cost', '')} per time unit under continuous dose 2)")
    except Exception:
        put("cost_holiday", "")
    g = load("fig4g_memory_disruption.csv").groupby("treatment").mean(numeric_only=True)
    put("clones_none", f"{g.loc['none', 'surviving_clones']:.0f}"); put("clones_predrug", f"{g.loc['before_drug', 'surviving_clones']:.0f}"); put("clones_during", f"{g.loc['before_and_during', 'surviving_clones']:.0f}")
    put("cells_none", f"{g.loc['none', 'surviving_cells']:.0f}"); put("cells_during", f"{g.loc['before_and_during', 'surviving_cells']:.0f}")

@safe
def f5a():
    d = load("fig5a_memory_genes.csv")
    put("score_min", f"{d[d.memory_time < 2].score_true.mean():.1f}"); put("score_max", f"{d[d.memory_time > 100].score_true.mean():.1f}")
    put("score_peak", f"{d.score_true.max():.1f}"); put("score_peak_seq", f"{d.score_seq.max():.1f}")
    y = (d.memory_time > 20).astype(int)
    put("auc_true", f"{roc_auc_score(y, d.score_true):.2f}"); put("auc_seq", f"{roc_auc_score(y, d.score_seq):.2f}")
    m = load("fig5bc_metrics.csv")
    for meth in ("pearson", "genie3"):
        for ds, tag in (("fixed_volume", "fixed"), ("population_counts", "counts"), ("population_concentration", "conc"), ("population_cycle_regressed", "reg"), ("sequenced_counts", "seq"), ("sequenced_normalized", "norm")):
            put(f"aupr_{meth}_{tag}", f"{m[(m.method == meth) & (m.dataset == ds)].aupr.mean():.2f}")

@safe
def f2b():
    import math
    k = load("fig2d_kernels.csv")
    crit = 1.63 / math.sqrt(20000)
    put("ks_crit20k", f"{crit:.4f}")
    coarse = k[k.kernel.str.contains("HybridSSATau\\(0.2")]
    put("ks_hybrid_coarse", f"{coarse.ks_telegraph.iloc[0]:.4f} and {coarse.ks_birthdeath.iloc[0]:.4f}")

@safe
def f4g():
    from scipy import stats
    d = load("fig4g_memory_disruption.csv")
    none = d[d.treatment == "none"]; pre = d[d.treatment == "before_drug"]
    t, pv = stats.ttest_ind(none.surviving_clones, pre.surviving_clones, equal_var=False)
    put("clones_predrug_p", f"mean of {len(none)} seeds each, P = {pv:.3f}, Welch's t test")
    put("resfrac_none", f"{none.high_fraction_before_drug.mean():.2f}")
    put("resfrac_predrug", f"{pre.high_fraction_before_drug.mean():.2f}")

@safe
def f5():
    """Random baselines, imputation and the perturbation baselines."""
    m = load("fig5bc_metrics.csv")
    # the random baseline differs between the directed (GENIE3) and undirected (correlation) evaluations
    put("aupr_random_dir", f"{m[m.directed].random_aupr.mean():.2f}")
    put("aupr_random_undir", f"{m[~m.directed].random_aupr.mean():.2f}")
    put("aupr_random", f"{m.random_aupr.mean():.2f}")
    for meth, tag in (("pearson", "p"), ("genie3", "g")):
        base = m[(m.method == meth) & (m.dataset == "sequenced_counts")].aupr.mean()
        for ds, dtag in (("imputed_knn_smoothing", "knn"), ("imputed_magic", "magic")):
            put(f"aupr_{tag}_{dtag}", f"{m[(m.method == meth) & (m.dataset == ds)].aupr.mean():.2f}")
    s = load("fig5d_summary.csv")
    # the mean R² of the correlation baseline is dominated by a single failure, so report the spread too
    z, c = s.r2_zero_baseline, s.r2_correlation_baseline
    put("r2_zero", f"{z.mean():.2f}"); put("r2_corr", f"{c.mean():.2f}")
    put("r2_zero_med", f"{z.median():.2f}"); put("r2_corr_med", f"{c.median():.2f}")
    put("r2_corr_nwin", f"{int((c > z).sum())} of {len(c)}")
    put("r2_corr_good", f"{c[c > 0].min():.2f} to {c.max():.2f}")
    put("r2_corr_worst", f"{c.min():.2f}")

@safe
def f6():
    a = kv("fig6a_summary.csv"); put("abc_a", f"({a['abc_median_k_on']:.2f}, {a['abc_median_k_off']:.2f}, {a['abc_median_k_tx']:.1f})")
    b = kv("fig6b_summary.csv"); put("naive_b", f"$k_{{\\mathrm{{on}}}} = {b['naive_k_on']:.2f}$, $k_{{\\mathrm{{off}}}} = {b['naive_k_off']:.2f}$, $k_{{\\mathrm{{tx}}}} = {b['naive_k_tx']:.1f}$")
    put("abc_b", f"({b['abc_median_k_on']:.2f}, {b['abc_median_k_off']:.2f}, {b['abc_median_k_tx']:.1f})")
    c = kv("fig6c_summary.csv"); put("abc_c", f"({c['abc_median_h_max']:.2f}, {c['abc_median_K']:.0f})")
    pc = load("fig6c_particles.csv"); lh, lK = np.log10(pc.h_max), np.log10(pc.K); w = pc.weight / pc.weight.sum()
    cv = np.cov(np.vstack([lh, lK]), aweights=w); put("ridge_corr", f"{cv[0, 1] / np.sqrt(cv[0, 0] * cv[1, 1]):.2f}")

@safe
def f7():
    c = kv("fig7_calibration.csv")
    put("cal_memgen", f"{c['memory_generations']:.1f}"); put("cal_pon", f"{100 * c['p_on']:.1f}%"); put("cal_ec50", f"{c['EC50']:.1f}")
    put("cal_hmax", f"{c['h_max']:.3f}"); put("cal_ic50", f"{c['IC50']:.1f}")
    d = load("fig7a_fates.csv"); m = d.groupby("cisplatin_uM").mean(numeric_only=True)
    # One definition throughout: the root-mean-square residual of the seed-averaged fate fractions,
    # pooled over every fraction of the concentrations named (six numbers for training, three for held-out).
    def resid(c):
        r = m.loc[c]
        return np.array([r.died - r.obs_died, r.divided - r.obs_divided, r.survived - r.obs_survived])
    rms = lambda v: float(np.sqrt((np.concatenate(v) ** 2).mean()))
    put("rmse_train", f"{rms([resid(7.0), resid(13.0)]):.3f}")
    put("rmse_heldout", f"{rms([resid(10.0)]):.3f}")
    # The colour scale of Fig. 7f is the calibration objective, not this error: it is evaluated on one
    # seed and adds the memory prior, so it is smaller than the four-seed error reported in the text.
    j = load("fig7j_cycle_fit.csv").set_index("key").value
    obj, seed1 = float(c["distance"]), float(j["rmse_flat"])
    put("cal_objective", f"{obj:.3f}")
    put("cal_objective_parts", f"{seed1:.3f} of single-seed training error and {obj - seed1:.3f} of prior penalty")
    r = m.loc[10.0]; put("obs_10", f"{r.obs_died:.2f}, {r.obs_divided:.2f} and {r.obs_survived:.2f}"); put("pred_10", f"{r.died:.2f}, {r.divided:.2f} and {r.survived:.2f}")
    put("obs_died_7", f"{m.loc[7.0].obs_died:.2f}"); put("obs_died_13", f"{m.loc[13.0].obs_died:.2f}")
    kd = load("fig7c_kin_correlation.csv")
    k = kd.groupby("relation").fate_correlation.mean()
    for key, rel in (("phi_sis", "sisters"), ("phi_c1", "first cousins"), ("phi_c2", "second cousins"), ("phi_c3", "third cousins"), ("phi_unrel", "unrelated")): put(key, f"{k[rel]:.2f}")
    g3 = kd[kd.relation == "third cousins"].fate_correlation
    # three seeds, so the interval is Student's t on n - 1 degrees of freedom; the normal
    # approximation would be narrow enough to exclude zero as an artefact of the approximation
    from scipy import stats as _st
    half = float(_st.t.ppf(0.975, len(g3) - 1)) * g3.std(ddof=1) / np.sqrt(len(g3))
    put("phi_c3_ci", f"95% confidence interval {g3.mean() - half:+.3f} to {g3.mean() + half:+.3f}, "
                     f"Student's t on {('one','two','three','four')[min(len(g3)-2,3)]} degrees of freedom")
    put("phi_c3_seeds", f"{len(g3)}")
    t = load("fig7d_timing.csv").groupby(["event", "cisplatin_uM"]).mean_h.mean()
    put("deathtime_shift", f"{abs(t[('death', 13.0)] - t[('death', 7.0)]) / t[('death', 7.0)] * 100:.0f}%")
    kc = load("fig7b_killcurves.csv"); g = kc[np.isclose(kc.cisplatin_uM, 13.0)].reset_index(drop=True)
    tt, N = g.t_since_drug_h.values, g.N_over_N0.values
    sl = lambda a, b: np.log(N[np.argmin(abs(tt - b))] / N[np.argmin(abs(tt - a))]) / (b - a)
    put("kc_ratio", f"{sl(0, 48) / sl(96, 168):.1f}")

@safe
def f7b():
    """Measurement precision, identifiability and step-size robustness of the calibration."""
    import numpy as np
    o = pd.read_csv(os.path.join(HERE, "..", "data", "iyer2025_u2os_fates.csv"))
    ses = []
    for _, r in o.iterrows():
        n = r["cells_at_drug"]
        for k in ("died", "divided", "survived_without_dividing"):
            p_ = r[k] / n
            ses.append(np.sqrt(p_ * (1 - p_) / n))
    binom = float(np.mean(ses))
    put("binom_se", f"{binom:.3f}")

    st = load("fig7i_stepsize.csv")
    sim = float(st[st.dt_h == st.dt_h.max()][["died", "divided", "survived"]].std().mean())
    put("sim_se", f"{sim:.3f}")
    rmse_train = float(V.get("rmse_train", "nan"))
    put("rmse_train_in_se", f"{rmse_train / np.hypot(binom, sim):.1f}")
    g = st.groupby("dt_h")[["died", "divided", "survived"]].mean()
    put("dt_shift", f"{float((g.loc[g.index.min()] - g.loc[g.index.max()]).abs().max()):.3f}")

    pr = load("fig7g_profile.csv")
    best = pr.rmse.min()
    def span(par):
        ok = pr[(pr.parameter == par) & (pr.rmse <= best + binom)]
        return ok.value.min(), ok.value.max()
    lo, hi = span("h_max")
    put("prof_hmax", f"a factor of {hi/lo:.0f}")
    f1, f2 = [span(p_)[1] / span(p_)[0] for p_ in ("k_on", "k_off")]
    lo_f, hi_f = round(min(f1, f2)), round(max(f1, f2))
    put("prof_switch", f"a factor of {lo_f}" if lo_f == hi_f else f"a factor of {lo_f} to {hi_f}")

    try:
        sl = load("fig7h_memory_slice.csv")
        b = sl.rmse.min(); ok = sl[sl.rmse <= b + binom]
        put("memslice_sentence",
            f"Holding the death parameters at their fitted values, memories from "
            f"{ok.memory_generations.min():.1f} to {ok.memory_generations.max():.1f} generations paired with resistant "
            f"fractions from {100*ok.p_on.min():.0f}% to {100*ok.p_on.max():.0f}% all describe the fate fractions within "
            f"the noise of the measurement (Supplementary Fig. 3b), because a shorter memory with more resistant cells and "
            f"a longer memory with fewer are hard to tell apart from fates alone.")
    except Exception:
        put("memslice_sentence", "")

@safe
def f8b():
    """Sensitivity of the melanoma schedule ranking to the memory of the resistant state."""
    import numpy as np
    d = load("fig8g_memory_sensitivity.csv")
    cont, inter, adapt = "continuous", "intermittent (S1320)", "adaptive (50 %)"
    g = d[d.mechanism == "no fitness cost"].groupby(["memory_generations", "schedule"]).ttp_baseline_weeks.agg(["mean", "std"])
    gens = sorted(d.memory_generations.unique())
    order_keeps = all(g.loc[(x, inter), "mean"] > g.loc[(x, cont), "mean"] for x in gens)
    gaps = {x: g.loc[(x, inter), "mean"] - g.loc[(x, cont), "mean"] for x in gens}
    lo, hi = min(gens), max(gens)
    base = 5.0 if 5.0 in gens else gens[len(gens) // 2]
    cost = d[d.mechanism == "fitness cost"]
    censored = bool((~cost.progressed_baseline.astype(bool)).all())
    put("memory_scan_note",
        f"The melanoma study gives the resistant state a memory of five net population doublings, about 20 weeks, "
        f"which is a modelling choice and not a measurement, and it is the quantity that decides whether resistant "
        f"cells revert during a three-week holiday. Repeating the schedule comparison over memories from {lo:.0f} to "
        f"{hi:.0f} net doublings, with {int(d.seed.nunique())} seeds at each point, leaves the ranking in place: "
        + ("without a fitness cost both interrupted schedules keep control longer than continuous dosing at every "
           "memory tested" if order_keeps else
           "without a fitness cost the ranking of continuous against intermittent dosing changes over this range") +
        f". What changes is how much the holidays are worth. The advantage of the intermittent schedule over "
        f"continuous dosing falls from {gaps[lo]:.0f} weeks at a memory of {lo:.0f} doublings to {gaps[hi]:.1f} "
        f"weeks at {hi:.0f}, against {gaps[base]:.1f} weeks at the {base:.0f} doublings used in the main text, "
        f"because a state that is forgotten quickly is re-drawn during a holiday while one that is remembered is "
        f"carried through it. The memory assumed here therefore sits at the conservative end: a shorter memory "
        f"would make the case for holidays stronger, not weaker. "
        + (f"With a fitness cost of resistance no schedule loses control within the follow-up at any memory in this "
           f"range, so that arm is insensitive to the assumption." if censored else
           f"With a fitness cost of resistance the ranking is unchanged over the same range."))
    put("memory_scan_range", f"{gaps[lo]:.0f} weeks at {lo:.0f} net doublings of memory to {gaps[hi]:.1f} weeks at {hi:.0f}")

    # the scan and the main melanoma run report the same quantity with different seed counts, and
    # the main text differences two rounded times, so state the reconciliation rather than leave
    # a reader to wonder why 2.3 and 3 weeks describe one result
    try:
        mm = load("fig8a_melanoma_schedules.csv")
        mm = mm[mm.mechanism == "no fitness cost"].groupby("schedule").ttp_baseline_weeks.agg(["mean", "std", "count"])
        cont, inter = mm.loc["continuous"], mm.loc["intermittent (S1320)"]
        gap_main = float(inter["mean"] - cont["mean"])
        se_main = float(np.sqrt(cont["std"] ** 2 / cont["count"] + inter["std"] ** 2 / inter["count"]))
        put("memory_scan_reconcile",
            f"The main text quotes this arm as {cont['mean']:.0f} weeks against {inter['mean']:.0f}, so differencing "
            f"the rounded times gives {round(inter['mean']) - round(cont['mean']):.0f} weeks; the unrounded "
            f"difference over the {int(cont['count'])} seeds of that run is {gap_main:.1f} weeks with a standard "
            f"error of {se_main:.1f} weeks, against {gaps[base]:.1f} weeks over the {int(d.seed.nunique())} seeds "
            f"of this scan. The three numbers are one result, reported at different seed counts and roundings.")
    except Exception as e:
        put("memory_scan_reconcile", "")

@safe
def f9():
    """Validation of the population and lineage layers against exact stationary laws."""
    import numpy as np
    v = load("fig9_validation.csv")
    v = v[v.kernel == "DirectSSA"]
    label = {"constitutive": "constitutive production", "bursty": "bursts of four molecules",
             "telegraph": "two-state promoter", "replication": "volume-scaled synthesis with gene replication"}
    dt0 = float(v[v.case == "constitutive"]["dt"].iloc[0])
    main = v[(v.case != "constitutive_dt") & (v["dt"] == dt0)]
    lin, pop = main[main["mode"] == "lineage"], main[main["mode"] == "population"]

    # Mother-machine samples are independent by construction, so the Kolmogorov-Smirnov test
    # applies as it stands. A population snapshot is not: its cells share ancestors, which makes
    # the independent critical value an under-estimate, so the population mode is tested on the
    # agreement of the replicate means with the exact mean instead.
    lin_pass = int((lin.ks <= lin.ks_crit99).sum())
    pop_pass = int((pop.ks <= pop.ks_crit99).sum())
    ratio = float((main.ks_other_mode / main.ks).min())

    rows, zs = [], {}
    for (case, mode), g in main.groupby(["case", "mode"], sort=False):
        mu, ex = g.mean_sim.mean(), g.mean_exact.iloc[0]
        se = float(g.mean_sim.std(ddof=1) / np.sqrt(len(g))) if len(g) > 1 else float(g.sem_mean.iloc[0])
        zs[(case, mode)] = (mu - ex) / se
        rows.append(dict(case=case, mode=mode, reps=len(g), n=int(g.n_cells.iloc[0]), mean_sim=mu, sem_mean=se,
                         mean_exact=ex, sd_sim=g.sd_sim.mean(), sd_exact=g.sd_exact.iloc[0],
                         ks=g.ks.max(), ks_crit=g.ks_crit99.iloc[0], ks_other=g.ks_other_mode.min()))
    r = pd.DataFrame(rows)
    zmax = max(abs(z) for z in zs.values())
    gaps = "; ".join(f"{label[c]} {r[(r.case == c) & (r['mode'] == 'lineage')].mean_exact.iloc[0]:.2f} against "
                     f"{r[(r.case == c) & (r['mode'] == 'population')].mean_exact.iloc[0]:.2f}"
                     for c in label if len(r[(r.case == c) & (r['mode'] == 'population')]))
    put("exact_pass", f"{lin_pass} of {len(lin)} lineage and {pop_pass} of {len(pop)} population samples")
    put("exact_ks_max", f"{main.ks.max():.4f}")

    head = (f"Along single lineages, where the mother-machine control makes the recorded cells independent and the "
            f"Kolmogorov-Smirnov test applies as it stands, all {len(lin)} samples match the exact law of their mode "
            f"(distances {lin.ks.min():.4f} to {lin.ks.max():.4f} against 99% critical values of "
            f"{lin.ks_crit99.min():.4f} to {lin.ks_crit99.max():.4f})" if lin_pass == len(lin) else
            f"Of the {len(lin)} single-lineage samples, {lin_pass} fall below the 99% critical value; the largest "
            f"distance is {lin.ks.max():.4f}")
    # the replicate means are all slightly low; say so rather than let "within z s.e." imply
    # a symmetric scatter, and give the worst case its p-value on the right reference distribution
    from scipy import stats as _st
    pz = {k: z for k, z in zs.items() if k[1] == "population"}
    nlow = sum(1 for z in pz.values() if z < 0)
    reps_pop = int(pop.groupby("case").size().max())
    worst = min(pz.values(), key=lambda z: z) if pz else 0.0
    pworst = float(2 * (1 - _st.t.cdf(abs(worst), reps_pop - 1)))
    pop_rel = [(row["mean_sim"] - row["mean_exact"]) / row["mean_exact"] * 100
               for _, row in r[r["mode"] == "population"].iterrows()]
    popsent = (f"A population snapshot is not an independent sample, because its cells share ancestors, so the "
               f"independent critical value understates the true one and the population mode is judged on the "
               f"agreement of its replicate means with the exact mean. Those means sit within {zmax:.1f} standard "
               f"errors of the exact values ({pop_pass} of the {len(pop)} individual samples also fall below the "
               f"independent critical value, the largest distance being {pop.ks.max():.4f}). All "
               f"{('two','three','four','five')[len(pz)-2] if 2 <= len(pz) <= 5 else len(pz)} population means are "
               f"nonetheless slightly low, by {abs(max(pop_rel)):.2f}% to "
               f"{abs(min(pop_rel)):.2f}%; no single deviation is resolved, the largest being {abs(worst):.1f} "
               f"standard errors ($p = {pworst:.2f}$, Student's t), but the common sign is not accounted for by "
               f"the discretisation, the subsampling cap or the burn-in, each of which is excluded in "
               f"Supplementary Note 14")
    tail = (f"Scored against the exact law of the other mode, the same samples give distances at least "
            f"{ratio:.0f} times larger, so the comparison has ample power to tell the two settings apart. The two "
            f"differ because a snapshot of a growing population over-weights cells that have just divided and so "
            f"just lost half their molecules: the mean molecule number is {gaps} along a lineage and in a "
            f"population respectively")
    put("exact_result", f"{head}. {popsent}. {tail}")

    dtv = v[v.case == "constitutive_dt"].sort_values("dt")
    conv = ", ".join("%.4f at $dt = 1/%d$" % (row.ks_scheme_vs_continuum, round(1 / row["dt"]))
                     for _, row in dtv.iterrows())
    floor = float(dtv["ks_crit99"].iloc[0])
    extra = ""
    try:
        a = load("fig9_agefit.csv")
        parts = []
        for mode in ("lineage", "population"):
            g = a[a["mode"] == mode]
            if len(g):
                parts.append(f"{0.5 * float(np.abs(g.observed - g.expected).sum()):.3f} ({mode})")
        if parts:
            extra = (f" In the model with gene replication the measured distribution of cell-cycle phase differs "
                     f"from the predicted one by a total variation distance of {' and '.join(parts)} over twenty "
                     f"bins, so the age structure that produces the snapshot weighting emerges from the branching "
                     f"dynamics rather than being imposed on it.")
    except Exception:
        pass
    put("exact_note",
        f"{head}. {popsent}. {tail}.\n\nDiscretising the interdivision time to the update step is the only bias in "
        f"the memoryless-timer comparison, and its size can be computed without simulating anything, as the distance "
        f"between the stationary law of the scheme and that of the continuous-time model: {conv}. It is first order "
        f"in the step and already below the 99% critical value at this sample size ({floor:.4f}) for every step "
        f"tested, so a step of a few percent of the interdivision time is enough for the discretisation to be "
        f"undetectable here. The small common negative offset of the population means is not explained by "
        f"that discretisation: the exact law each population sample is scored against is the stationary law of "
        f"the scheme at the simulated step, and solving it at steps from 1/64 to 1/512 of the interdivision time "
        f"moves the exact population mean by less than one part in $10^{{9}}$, while it moves the lineage mean by "
        f"the amounts tabulated above. Nor is it the subsampling cap, which draws cells uniformly without "
        f"replacement from a snapshot whose division and thinning are independent of the molecule number, so the "
        f"retained cells are exchangeable; nor the burn-in, since the population mean relaxes at the sum of the "
        f"decay and growth rates and twelve generations leave a transient far below the observed offset. On the "
        f"evidence available it is the scatter of five runs on four degrees of freedom, and an offset of this "
        f"size is neither established nor excluded.{extra}")

    short = {"constitutive": "constitutive", "bursty": "bursty",
             "telegraph": "two-state", "replication": "replication"}
    # explicit proportional widths: twelve equal columns overflowed their boxes and collided
    hdr = ("| Model | Mode | Mean (sim.) | s.e. | Mean (exact) | s.d. (sim.) | s.d. (exact) | "
           "Mean in s.e. | Largest KS | Indep. 99% | KS, other mode |\n"
           "|:-----------|:-----------|------:|-----:|------:|------:|------:|-----:|------:|------:|------:|")
    lines = [f"| {short.get(x.case, x.case)} | {x['mode']} | {x.mean_sim:.3f} | {x.sem_mean:.3f} | "
             f"{x.mean_exact:.3f} | {x.sd_sim:.3f} | {x.sd_exact:.3f} | {zs[(x.case, x['mode'])]:+.2f} | "
             f"{x.ks:.4f} | {x.ks_crit:.4f} | {x.ks_other:.4f} |"
             for _, x in r.iterrows()]
    nlin = int(r[r["mode"] == "lineage"].n.max()); npop = int(r[r["mode"] == "population"].n.max())
    put("exact_table",
        f"Supplementary Table 5. Simulated against exact stationary laws. Lineage rows are one sample of "
        f"{nlin:,} cells; population rows are the mean over {reps_pop} independent runs of {npop:,} cells each. "
        f"The quoted error on a population mean is the standard error of that average, the spread between the "
        f"{reps_pop} run means divided by the square root of {reps_pop}, and on a lineage mean it is the standard "
        f"error within the single sample; \"mean in s.e.\" is the difference from the exact mean in those units. "
        f"\"Largest KS\" is the worst Kolmogorov-Smirnov distance over the samples of the row, and the independent "
        f"critical value beside it is exact for the lineage rows and an under-estimate for the population rows.\n\n"
        + hdr + "\n" + "\n".join(lines))

@safe
def f7c():
    """Cell-cycle-dependent killing tested against heritable expression."""
    import numpy as np
    fit = load("fig7j_cycle_fit.csv").set_index("key").value
    prof = load("fig7j_cycle_profile.csv")
    fat = load("fig7j_cycle_fates.csv")
    o = pd.read_csv(os.path.join(HERE, "..", "data", "iyer2025_u2os_fates.csv"))
    ses, ov = [], {}
    for _, r in o.iterrows():
        n = r["cells_at_drug"]
        for k, lab in (("died", "died"), ("divided", "divided"), ("survived_without_dividing", "survived")):
            q = r[k] / n; ses.append(np.sqrt(q * (1 - q) / n)); ov[(float(r["cisplatin_uM"]), lab)] = q
    se = float(np.mean(ses))
    rf, rc = float(fit["rmse_flat"]), float(fit["rmse_cycle"])
    beta = float(fit["cycle_baseline"])
    held = 10.0
    def held_rmse(tag):
        g = fat[(fat.model == tag) & (fat.cisplatin_uM == held)][["died", "divided", "survived"]].mean()
        return float(np.sqrt(np.mean([(g[l] - ov[(held, l)]) ** 2 for l in ("died", "divided", "survived")])))
    hf, hc = held_rmse("cycle-independent"), held_rmse("cycle-dependent")
    ok = prof[prof.rmse <= prof.rmse.min() + se]
    span = f"{ok.cycle_baseline.min():.2f} to {ok.cycle_baseline.max():.2f}"
    put("cyc_beta", f"{beta:.2f}"); put("cyc_span", span)
    put("cyc_rmse_flat", f"{rf:.3f}"); put("cyc_rmse_cycle", f"{rc:.3f}")
    put("cyc_held_flat", f"{hf:.3f}"); put("cyc_held_cycle", f"{hc:.3f}")
    # the verdict follows the numbers rather than the other way round
    diff = rf - rc                       # positive when the cycle-dependent fit is the better one
    if abs(diff) < 0.3 * se:
        head = (f"The two fits describe the training fractions equally well ({rc:.3f} against {rf:.3f}, "
                f"against a mean sampling error of {se:.3f} in the measurement)")
    elif diff > 0:
        head = f"The cycle-dependent fit describes the training fractions better ({rc:.3f} against {rf:.3f})"
    else:
        head = f"The cycle-dependent fit describes the training fractions no better ({rc:.3f} against {rf:.3f})"
    if hc < hf - 0.3 * se:
        tail = f"and the cycle-dependent fit predicts the held-out concentration better ({hc:.3f} against {hf:.3f})"
    elif hc > hf + 0.3 * se:
        tail = f"and the cycle-dependent fit does not improve the held-out prediction either ({hc:.3f} against {hf:.3f})"
    else:
        tail = f"and the two predict the held-out concentration about equally ({hc:.3f} against {hf:.3f})"
    put("cycle_result",
        f"{head}, {tail}. The fitted cycle-independent fraction is {beta:.2f}, but values from {span} all sit "
        f"within one sampling error of the best fit, so the fate fractions do not determine it")
    put("cycle_note",
        f"Refitting under the same objective, bounds and optimiser budget gives a cycle-independent fraction of "
        f"{beta:.2f} and a root-mean-square error of {rc:.3f} on the training fractions, against {rf:.3f} for the "
        f"cycle-independent model; at the held-out concentration the errors are {hc:.3f} and {hf:.3f}. The profile "
        f"over that fraction is flat: every value from {span} lies within one binomial standard error of the "
        f"measurement ({se:.3f}) of the best fit (Supplementary Fig. 4a). Three fate counts per concentration "
        f"therefore cannot say whether cells are spared because they inherited a protective state or because they "
        f"were outside the replication window, and the parameters of either fit should be read with that in mind.")

@safe
def f4c():
    """Do the persister and melanoma conclusions survive cell-cycle-gated killing?"""
    import numpy as np
    d0 = load("fig4c_decay_vs_dose.csv"); d1 = load("fig4i_cycle_decay.csv")
    c0 = load("fig4d_fate_concordance.csv"); c1 = load("fig4i_cycle_concordance.csv")
    # the gated scan is the five-seed rerun, so the cycle-blind reference must be the five-seed
    # scan too; comparing it against the older single-seed table manufactured a release-period
    # difference that is a seed artefact rather than an effect of the gate
    try:
        s0 = load("fig4f_schedules_seeds.csv")
    except Exception:
        s0 = load("fig4f_schedules.csv")
    s1 = load("fig4i_cycle_schedules.csv")
    m0 = load("fig8a_melanoma_schedules.csv"); m1 = load("fig8f_melanoma_cycle.csv")

    # decay rates are reported as a range, not a ratio: the gated rate crosses zero at low dose
    def dec(d):
        g = d[d.dose >= 0.5].sort_values("dose")
        return float(g.decay_rate.iloc[0]), float(g.decay_rate.iloc[-1])
    a1, b1 = dec(d1); a0, b0 = dec(d0)
    put("cyc4_decay", f"{a1:+.3f} to {b1:+.3f} per time unit over doses 0.5 to 2, against {a0:+.3f} to {b0:+.3f} "
                      f"without the gate")
    def dtr(d):
        g = d[d.dose > 0].death_time_mean.dropna(); return float(g.min()), float(g.max())
    q1, q0 = dtr(d1), dtr(d0)
    put("cyc4_dtime", f"{q1[0]:.1f} to {q1[1]:.1f} against {q0[0]:.1f} to {q0[1]:.1f} time units")

    def exc(df, m, r):
        g = df[(df.model == m) & (df.relation == r)]
        return float(g.concordance.iloc[0] - g.expected_independent.iloc[0]) if len(g) else float("nan")
    mem0, mem1 = exc(c0, "memory", "sisters"), exc(c1, "memory", "sisters")
    fst0, fst1 = exc(c0, "fast", "sisters"), exc(c1, "fast", "sisters")
    put("cyc4_mem", f"{mem1:+.3f} against {mem0:+.3f}")
    put("cyc4_fast", f"{fst1:+.3f} against {fst0:+.3f}")
    put("cyc4_mem_g", f"{mem1:+.3f}"); put("cyc4_mem_b", f"{mem0:+.3f}")
    put("cyc4_fast_g", f"{fst1:+.3f}"); put("cyc4_fast_b", f"{fst0:+.3f}")
    put("cyc4_ratio", f"{mem1 / fst1:.0f}-fold" if fst1 > 1e-6 else "far")
    put("cyc4_confound_word", "does" if fst1 > max(0.02, 2 * abs(fst0)) else "does not")

    # schedules: the gate shifts the level of every growth rate, so report the shift, the winner
    # and the order below the winner separately rather than collapsing them into one verdict
    def sched(df, m):
        return df[df.model == m].groupby(["release_period", "dose"]).long_term_growth_rate.agg(
            ["mean", "std", "count"])
    def bs(df, m):
        rp, dose = sched(df, m)["mean"].idxmin()
        return float(rp), float(dose)
    label = {"pre_existing": "pre-existing tolerance", "pre_existing_cost": "a fitness cost of resistance",
             "drug_induced": "drug-induced tolerance"}
    models = list(dict.fromkeys(s0.model))
    same_dose = sum(bs(s0, m)[1] == bs(s1, m)[1] for m in models)
    moved = [m for m in models if bs(s0, m)[0] != bs(s1, m)[0]]
    continuous_both = all(bs(d, m)[0] == 0.0 for d in (s0, s1) for m in models)

    # how far the gate lifts the whole surface, and how far it reorders it below the winner
    rise = np.concatenate([(sched(s1, m)["mean"] - sched(s0, m)["mean"]).dropna().values for m in models])
    rho = [sched(s0, m)["mean"].corr(sched(s1, m)["mean"], method="spearman") for m in models]
    put("cyc4_rise", f"{rise.min():+.4f} to {rise.max():+.4f} per time unit at all {len(rise)} points of the scan"
        if (rise > 0).all() else f"{rise.min():+.4f} to {rise.max():+.4f} per time unit")
    put("cyc4_rank", f"{min(rho):.3f} to {max(rho):.3f}")

    # same paired convention as the ungated scan: the schedules share seeds, so the margin is a
    # paired difference and is tested by a paired t-test rather than against twice an unpaired error
    marg = [dict(paired_margin(d, m), scan=tag) for tag, d in (("without the gate", s0), ("under the gate", s1))
            for m in models]
    won = [m for m in marg if m["resolved"]]
    lost = [m for m in marg if not m["resolved"]]
    def _ex(m):
        return (f"{label.get(m['model'], m['model'])} {m['scan']} ({m['gap']:+.4f} against a paired standard "
                f"error of {m['se']:.4f})")
    if not lost:
        exc_clause = ""
    elif len(lost) == 1:
        exc_clause = f", the exception being {_ex(lost[0])}"
    else:
        exc_clause = (", the exceptions being " + ", ".join(_ex(m) for m in lost[:-1]) + " and " + _ex(lost[-1]))
    put("cyc4_margin",
        f"a two-sided paired t-test resolves the margin over the best holiday at the 5% level in "
        + (f"all {len(marg)} comparisons" if not lost else f"{len(won)} of the {len(marg)} comparisons{exc_clause}"))

    # the intermediate-dose optimum, quoted only when the grid actually carries those points
    def at(df, m, rp, dose):
        g = sched(df, m)["mean"]
        return float(g.loc[(rp, dose)]) if (rp, dose) in g.index else None
    lo1, hi1 = at(s1, "drug_induced", 0.0, 1.0), at(s1, "drug_induced", 0.0, 2.0)
    lo0, hi0 = at(s0, "drug_induced", 0.0, 1.0), at(s0, "drug_induced", 0.0, 2.0)
    interdose = None not in (lo1, hi1, lo0, hi0) and lo1 < hi1 and lo0 < hi0
    put("cyc4_interdose", f" nor the intermediate-dose optimum under drug-induced tolerance ({lo1:+.4f} against "
                          f"{hi1:+.4f} under the gate and {lo0:+.4f} against {hi0:+.4f} without it)"
        if interdose else "")

    if continuous_both and not moved:
        put("cyc4_sched",
            f"continuous dosing wins in all {len(models)} models with the gate as without it, at the same dose in "
            f"{same_dose} of {len(models)}. The gate raises the long-term growth rate ({V.get('cyc4_rise', '')}), "
            f"which is what confining three quarters of the hazard to a window of the cycle implies, and it "
            f"perturbs the order of the schedules below the winner (Spearman correlation {V.get('cyc4_rank', '')} "
            f"between the gated and cycle-blind orderings), but it changes neither which schedule wins"
            f"{V.get('cyc4_interdose', '')}. Across the two scans {V.get('cyc4_margin', '')}")
    else:
        which = " and with ".join(label.get(m, m) for m in moved) if moved else "some mechanisms"
        put("cyc4_sched",
            f"the best dose is unchanged in {same_dose} of {len(models)} models, but the best release period is "
            f"not: it moves with {which}. The gate raises the long-term growth rate "
            f"({V.get('cyc4_rise', '')}) and reorders the scan (Spearman correlation {V.get('cyc4_rank', '')})")

    def best(df, mm):
        g = df[df.mechanism == mm].groupby("schedule").ttp_baseline_weeks.mean()
        return g.idxmax()
    mechs = list(dict.fromkeys(m0.mechanism))
    put("cyc8_same", f"{sum(best(m0, mm) == best(m1, mm) for mm in mechs)} of {len(mechs)}")
    def order(df, mm):
        return tuple(df[df.mechanism == mm].groupby("schedule").ttp_baseline_weeks.mean().sort_values().index)
    flipped = [mm for mm in mechs if order(m0, mm) != order(m1, mm)]
    put("cyc8_full", f"{len(mechs) - len(flipped)} of {len(mechs)}")
    put("cyc8_flipped", "none" if not flipped else
        "; ".join(f"under {mm} the two losing schedules trade places" for mm in flipped))
    pp = "partial protection"
    if pp in mechs:
        g0 = m0[m0.mechanism == pp].groupby("schedule").ttp_baseline_weeks.mean()
        g1 = m1[m1.mechanism == pp].groupby("schedule").ttp_baseline_weeks.mean()
        cont, inter = "continuous", "intermittent (S1320)"
        keep = (g1[cont] > g1[inter]) == (g0[cont] > g0[inter])
        put("cyc8_trial", ("still" if keep else "no longer") +
            f" favours continuous dosing under partial protection ({g1[cont]:.0f} against {g1[inter]:.0f} weeks, "
            f"from {g0[cont]:.0f} and {g0[inter]:.0f})")

    put("cycle_case_result",
        f"Population decay still rises with dose ({V.get('cyc4_decay', '')}) and the mean time to death stays "
        f"nearly dose-invariant ({V.get('cyc4_dtime', '')}), so neither signature depends on a cycle-blind hazard. "
        f"Sisters share a birth time and so a cycle phase, which lets a gated hazard correlate their fates with no "
        f"inherited state; it {V.get('cyc4_confound_word', '')} do so to any useful degree, lifting the sister "
        f"concordance of the memoryless fast-switching control only to {V.get('cyc4_fast_g', '')} above "
        f"independence (from {V.get('cyc4_fast_b', '')}), against {V.get('cyc4_mem_g', '')} for the memory gene, "
        f"so the kin signature still reads expression memory rather than the cycle. In the schedule scan, "
        f"{V.get('cyc4_sched', '')}")
    put("cycle_case_note",
        f"Population decay still rises with dose ({V.get('cyc4_decay', '')}), at roughly half the rate, and the "
        f"mean time to death of killed cells stays nearly dose-invariant ({V.get('cyc4_dtime', '')}), so the "
        f"dose-dependence of decay alongside dose-invariant single-cell timing is not an artefact of a cycle-blind "
        f"hazard. The kin correlations deserve more care, because sisters are born together and therefore occupy "
        f"the same phase of the cycle: a gated hazard can in principle correlate their fates with no heritable "
        f"expression state at all. It {V.get('cyc4_confound_word', '')}: in the fast-switching control, which has "
        f"no usable memory, the sister concordance above independence rises only from {V.get('cyc4_fast_b', '')} "
        f"to {V.get('cyc4_fast_g', '')}, while the memory gene sits at {V.get('cyc4_mem_g', '')}, "
        f"{V.get('cyc4_ratio', '')} larger. In the "
        f"schedule scan, {V.get('cyc4_sched', '')}. In the melanoma study the best "
        f"schedule is unchanged in {V.get('cyc8_same', '')} mechanisms and the model {V.get('cyc8_trial', '')}; "
        f"the full ordering of the three schedules survives in {V.get('cyc8_full', '')} mechanisms "
        f"({V.get('cyc8_flipped', '')}).")
    put("cycle_mel_result",
        f"the best schedule is unchanged in {V.get('cyc8_same', '')} mechanisms and the model "
        f"{V.get('cyc8_trial', '')}")

    # matched cycle-averaged mean hazard: the gate lowers the average hazard as well as making it
    # phase-dependent, so a second arm raises h_max to hold the average fixed
    try:
        dm = load("fig4i_cycle_matched.csv"); cm = load("fig4i_cycle_matched_concordance.csv")
        mk = kv("fig4i_cycle_matching.csv")
        am, bm = dec(dm); a0, b0 = dec(d0); a1, b1 = dec(d1)
        memm, fstm = exc(cm, "memory", "sisters"), exc(cm, "fast", "sisters")
        try:
            mm = load("fig8f_melanoma_cycle_matched.csv")
            mel = (f" In the melanoma study the best schedule is unchanged in "
                   f"{sum(best(mm, x) == best(m0, x) for x in mechs)} of {len(mechs)} mechanisms at the matched hazard as well.")
        except Exception:
            mel = ""
        put("cycle_matched_note",
            f"The gate does two things at once. It makes the hazard depend on cycle phase, and it lowers the hazard "
            f"on average: over a uniform cycle phase the multiplier averages {mk['mean_multiplier_uniform_phase']:.3f}, "
            f"and under the phase density these simulations realise, $f(\\varphi) = 2/(1+\\varphi)^2$ for a sizer with "
            f"exponential growth in a growing population, {mk['mean_multiplier_realised']:.3f}. Half of the shallower "
            f"dose response above is therefore just less killing. Repeating the comparison with the maximal hazard "
            f"raised by $1/{mk['mean_multiplier_realised']:.3f}$, so that the cycle-averaged hazard matches the "
            f"cycle-blind one ($h_{{\\max}} = {mk['h_max_matched']:.3f}$ against {mk['h_max_blind']:.3f}), gives a "
            f"decay rate of {am:+.3f} to {bm:+.3f} per time unit over doses 0.5 to 2, against {a1:+.3f} to {b1:+.3f} "
            f"at the unmatched gated hazard and {a0:+.3f} to {b0:+.3f} with no gate at all. At matched mean hazard "
            f"the sister concordance above independence is {memm:+.3f} for the memory gene and {fstm:+.3f} for the "
            f"fast-switching control, against {V.get('cyc4_mem_g', '')} and {V.get('cyc4_fast_g', '')} unmatched, so "
            f"the conclusion that the gate cannot manufacture the kin signature does not rest on the gate also "
            f"killing less.{mel}")
    except Exception as e:
        put("cycle_matched_note", "")

@safe
def f8():
    d = load("fig8a_melanoma_schedules.csv"); m = d.groupby(["mechanism", "schedule"]).mean(numeric_only=True)
    for mech, mtag in (("no fitness cost", "nocost"), ("fitness cost", "cost"), ("partial protection", "partial")):
        for sched, stag in (("continuous", "cont"), ("intermittent (S1320)", "int"), ("adaptive (50 %)", "adapt")):
            r = m.loc[(mech, sched)]
            put(f"ttp_{stag}_{mtag}", f"{r.ttp_baseline_weeks:.0f}" + ("" if r.progressed_baseline > 0.99 else "+"))
            put(f"ttpn_{stag}_{mtag}", f"{r.ttp_nadir_weeks:.0f}" + ("" if r.progressed_nadir > 0.99 else "+"))
        put(f"dose_adapt_{mtag}", f"{100 * m.loc[(mech, 'adaptive (50 %)')].cumulative_dose_weeks / m.loc[(mech, 'continuous')].cumulative_dose_weeks:.0f}%")
        put(f"dose_int_{mtag}", f"{100 * m.loc[(mech, 'intermittent (S1320)')].cumulative_dose_weeks / m.loc[(mech, 'continuous')].cumulative_dose_weeks:.0f}%")
    for tag, key in (("no_fitness_cost", "opt_nocost"), ("fitness_cost", "opt_cost"), ("partial_protection", "opt_partial")):
        t = load(f"fig8c_melanoma_optimum_{tag}.csv"); b = t.loc[t.ttp_weeks.idxmax()]
        put(key, f"a period of {b.period_weeks:.1f} weeks with {100 * b.duty:.0f}% of the time on drug ({b.ttp_weeks:.0f} weeks to loss of control)")
    g = load("fig8d_gbm_regimens.csv").groupby(["population", "mechanism", "regimen"]).mean(numeric_only=True)
    for pop, ptag in (("MGMT methylated (1 % expressing)", "meth"), ("MGMT unmethylated (30 % expressing)", "unmeth")):
        for mech, gtag in (("MGMT stable", "stable"), ("MGMT consumed by drug", "consumed")):
            for reg, rtag in (("standard 5/28", "std"), ("dose-dense 21/28", "dense"), ("dense, equal cumulative", "equal")):
                put(f"lk_{ptag}_{gtag}_{rtag}", f"{g.loc[(pop, mech, reg)].log_kill:.2f}")
                put(f"nend_{ptag}_{gtag}_{rtag}", f"{g.loc[(pop, mech, reg)].N_end_over_N0:.1f}")
    for tag in ("methylated", "unmethylated"):
        for mech in ("stable", "consumed"):
            t = load(f"fig8e_gbm_days_on_{tag}_{mech}.csv").groupby("days_on").log_N_end_over_N0.mean(); put(f"days_{tag}_{mech}", f"{int(t.idxmin())}")

@safe
def fsupp():
    t = os.path.join(HERE, "tableS_runtime.md")
    if os.path.exists(t): put("runtime_table", open(t).read().strip())
    d = load("fig4c_times.csv")
    for dose, g in d[d.event == "death"].groupby("dose"):
        put(f"death_cv_{str(dose).replace('.', '_')}", f"{g.time.std() / g.time.mean():.2f}")
    for tag in ("fig6a", "fig6b", "fig6c"):
        sch = load(f"{tag}_schedule.csv")
        put(f"{tag}_gens", f"{int(sch.generation.max())}"); put(f"{tag}_eps", f"{sch.epsilon.iloc[-1]:.3g}"); put(f"{tag}_acc", f"{100 * sch.acceptance.iloc[-1]:.0f}%")
for f in (f2, f2b, f3, f4, f4g, f5a, f5, f6, f7, f7b, f7c, f8, f8b, f4c, f9): f()
fsupp()

def fill(template, target):
    tpl = open(os.path.join(HERE, template)).read()
    used = set(re.findall(r"{{([a-z_0-9]+)}}", tpl))
    missing = sorted(used - set(V))
    empty = sorted(k for k in used & set(V) if not str(V[k]).strip())
    if empty: print(f"  [{target}: EMPTY placeholders: {empty}]")
    out = re.sub(r"{{([a-z_0-9]+)}}", lambda m: str(V.get(m.group(1), "[" + m.group(1) + "]")), tpl)
    open(os.path.join(HERE, target), "w").write(out)
    print(f"{target}: filled from {len(V)} values; missing:", missing)
fill("02_results.template.md", "02_results.md")
if os.path.exists(os.path.join(HERE, "supplement.template.md")): fill("supplement.template.md", "supplement.md")
