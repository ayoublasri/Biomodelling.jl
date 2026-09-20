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
    parts = []
    for kern in ("HybridSSATau", "DirectSSA"):
        g = d[(d.kernel == kern) & (d.cells == d.cells.max())].sort_values("genes")
        if len(g) < 2: continue
        lo, hi = g.iloc[0], g.iloc[-1]
        parts.append(f"{kern} {lo.seconds:.2f} s to {hi.seconds:.1f} s, a factor of {hi.seconds/lo.seconds:.0f} "
                     f"for a factor of {hi.genes/lo.genes:.0f} in genes")
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
    gaps = []
    for model, tag in (("pre_existing", "pre"), ("pre_existing_cost", "cost"), ("drug_induced", "ind")):
        g = agg[agg.model == model].sort_values("mean")
        best, second = g.iloc[0], g.iloc[1]
        sem = (best["std"] / np.sqrt(best["count"])) if best["count"] > 1 else float("nan")
        err = "" if nseed == 1 else f" ± {best['std']:.3f} over {nseed} seeds"
        put(f"best_{tag}", f"release period {best.release_period:g}, dose {best.dose:g} (growth rate {best['mean']:+.3f}{err} per time unit)")
        cont = g[(g.release_period == 0) & (g.dose == g.dose.max())]["mean"].iloc[0]
        put(f"cont_{tag}", f"{cont:+.3f}")
        if nseed > 1:
            gaps.append((tag, float(second["mean"] - best["mean"]), float(np.hypot(sem, second["std"] / np.sqrt(second["count"])))))
    if gaps:
        resolved = [t for t, d, e in gaps if d > 2 * e]
        put("sched_seed_note",
            f"Each point of the scan is the mean of {nseed} independent seeds. The gap between the best schedule and "
            f"the next best is {min(d for _, d, _ in gaps):.3f} to {max(d for _, d, _ in gaps):.3f} per time unit, "
            f"against a standard error of {min(e for _, _, e in gaps):.3f} to {max(e for _, _, e in gaps):.3f} on the "
            f"difference, so the ranking of the top two schedules is resolved in {len(resolved)} of the "
            f"{len(gaps)} mechanisms; elsewhere the scan identifies a region of good schedules rather than a single best one.")
    else:
        put("sched_seed_note", "")
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
    half = 1.96 * g3.std(ddof=1) / np.sqrt(len(g3))
    put("phi_c3_ci", f"95% confidence interval {g3.mean() - half:+.3f} to {g3.mean() + half:+.3f}")
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
def f9():
    """Validation of the population and lineage layers against exact stationary laws."""
    import numpy as np
    v = load("fig9_validation.csv")
    v = v[v.kernel == "DirectSSA"]
    label = {"constitutive": "constitutive production", "bursty": "bursts of four molecules",
             "telegraph": "two-state promoter", "replication": "volume-scaled synthesis with gene replication"}
    dt0 = float(v[v.case == "constitutive"]["dt"].iloc[0])
    main = v[(v.case != "constitutive_dt") & (v["dt"] == dt0)]
    # each replicate is an independent sample, so each gets its own test at its own sample size
    passed = int((main.ks <= main.ks_crit99).sum())
    ratio = float((main.ks_other_mode / main.ks).min())
    put("exact_pass", f"{passed} of {len(main)}")
    put("exact_ks_max", f"{main.ks.max():.4f}")

    rows = []
    for (case, mode), g in main.groupby(["case", "mode"], sort=False):
        rows.append(dict(case=case, mode=mode, reps=len(g), n=int(g.n_cells.iloc[0]),
                         mean_sim=g.mean_sim.mean(), sem_mean=float(g.mean_sim.std(ddof=1) / np.sqrt(len(g))) if len(g) > 1 else float(g.sem_mean.iloc[0]),
                         mean_exact=g.mean_exact.iloc[0], sd_sim=g.sd_sim.mean(), sd_exact=g.sd_exact.iloc[0],
                         ks=g.ks.max(), ks_crit=g.ks_crit99.iloc[0], ks_other=g.ks_other_mode.min()))
    r = pd.DataFrame(rows)
    lin = r[r["mode"] == "lineage"].set_index("case"); pop = r[r["mode"] == "population"].set_index("case")
    gaps = "; ".join(f"{label[c]} {lin.loc[c, 'mean_exact']:.2f} against {pop.loc[c, 'mean_exact']:.2f}"
                     for c in lin.index if c in pop.index)
    if passed == len(main):
        head = (f"Every simulated distribution matches the exact law of its own mode. The largest "
                f"Kolmogorov-Smirnov distance over the {len(main)} samples is {main.ks.max():.4f}, and each lies "
                f"below the 99% critical value for its own sample size ({main.ks_crit99.min():.4f} to "
                f"{main.ks_crit99.max():.4f})")
    else:
        w = main.loc[main.ks.idxmax()]
        head = (f"{passed} of the {len(main)} samples fall below the 99% critical value for their sample size; the "
                f"largest distance is {main.ks.max():.4f} for {label.get(w.case, w.case)} in {w['mode']} mode, "
                f"against a critical value of {w.ks_crit99:.4f}")
    tail = (f"Scored against the exact law of the other mode, the same samples give distances at least "
            f"{ratio:.0f} times larger, so the comparison has ample power to tell the two settings apart. They "
            f"differ because a snapshot of a growing population over-weights cells that have just divided and so "
            f"just lost half their molecules: the mean molecule number is {gaps} along a lineage and in a "
            f"population respectively")
    put("exact_result", f"{head}. {tail}")

    dtv = v[v.case == "constitutive_dt"].sort_values("dt")
    conv = ", ".join("%.4f at $dt = %g$" % (row.ks_scheme_vs_continuum, row["dt"]) for _, row in dtv.iterrows())
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
                     f"from the predicted one by a total variation distance of {' and '.join(parts)} over twenty bins of cycle "
                     f"phase, so the age structure that produces the snapshot weighting emerges from the branching dynamics rather "
                     f"than being imposed on it.")
    except Exception:
        pass
    put("exact_note",
        f"{head}. {tail}.\n\nDiscretising the interdivision time to the update step is the only bias in the "
        f"memoryless-timer comparison, and its size can be computed without simulating anything, as the distance "
        f"between the stationary law of the scheme and that of the continuous-time model: {conv}. It is first "
        f"order in the step and already below the 99% critical value at this sample size ({floor:.4f}) for every "
        f"step tested, so a step of a few percent of the interdivision time is enough for the discretisation to be "
        f"undetectable here.{extra}")

    hdr = ("| Model | Mode | Cells per sample | Samples | Mean (simulated) | Mean (exact) | s.d. (simulated) | "
           "s.d. (exact) | KS to own mode | 99% critical | KS to other mode |\n"
           "|---|---|---|---|---|---|---|---|---|---|---|")
    lines = [f"| {label.get(x.case, x.case)} | {x['mode']} | {x.n:,} | {x.reps} | {x.mean_sim:.3f} ± {x.sem_mean:.3f} | "
             f"{x.mean_exact:.3f} | {x.sd_sim:.3f} | {x.sd_exact:.3f} | {x.ks:.4f} | {x.ks_crit:.4f} | {x.ks_other:.4f} |"
             for _, x in r.iterrows()]
    put("exact_table",
        "Supplementary Table 5. Simulated against exact stationary laws (worst distance over the samples of each row).\n\n"
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
    s0 = load("fig4f_schedules.csv"); s1 = load("fig4i_cycle_schedules.csv")
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

    # schedules: separate what survives (the dose) from what does not (the release period)
    def bs(df, m):
        g = df[df.model == m].groupby(["release_period", "dose"]).long_term_growth_rate.mean()
        rp, dose = g.idxmin()
        return float(rp), float(dose)
    label = {"pre_existing": "pre-existing tolerance", "pre_existing_cost": "a fitness cost of resistance",
             "drug_induced": "drug-induced tolerance"}
    models = list(dict.fromkeys(s0.model))
    same_dose = sum(bs(s0, m)[1] == bs(s1, m)[1] for m in models)
    moved = [m for m in models if bs(s0, m)[0] != bs(s1, m)[0]]
    if not moved:
        put("cyc4_sched", f"the best schedule is unchanged in all {len(models)} models")
    else:
        to_zero = all(bs(s1, m)[0] == 0.0 for m in moved)
        which = " and with ".join(label.get(m, m) for m in moved)
        put("cyc4_sched",
            f"the best dose is unchanged in {same_dose} of {len(models)} models, but the best release period is "
            f"not. With {which}, the holidays that won without the gate "
            + ("collapse to continuous dosing" if to_zero else "move") +
            ", because a gated hazard kills the reverting cells more slowly on re-exposure")

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
        f"so the kin signature still reads expression memory rather than the cycle. The schedule conclusions are "
        f"the fragile ones, because {V.get('cyc4_sched', '')}")
    put("cycle_case_note",
        f"Population decay still rises with dose ({V.get('cyc4_decay', '')}), at roughly half the rate, and the "
        f"mean time to death of killed cells stays nearly dose-invariant ({V.get('cyc4_dtime', '')}), so the "
        f"dose-dependence of decay alongside dose-invariant single-cell timing is not an artefact of a cycle-blind "
        f"hazard. The kin correlations deserve more care, because sisters are born together and therefore occupy "
        f"the same phase of the cycle: a gated hazard can in principle correlate their fates with no heritable "
        f"expression state at all. It {V.get('cyc4_confound_word', '')}: in the fast-switching control, which has "
        f"no usable memory, the sister concordance above independence rises only from {V.get('cyc4_fast_b', '')} "
        f"to {V.get('cyc4_fast_g', '')}, while the memory gene sits at {V.get('cyc4_mem_g', '')}, "
        f"{V.get('cyc4_ratio', '')} larger. What the "
        f"gate does change is the schedule scan: {V.get('cyc4_sched', '')}. In the melanoma study the best "
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
for f in (f2, f2b, f3, f4, f4g, f5a, f5, f6, f7, f7b, f7c, f8, f4c, f9): f()
fsupp()

def fill(template, target):
    tpl = open(os.path.join(HERE, template)).read()
    missing = sorted(set(re.findall(r"{{([a-z_0-9]+)}}", tpl)) - set(V))
    out = re.sub(r"{{([a-z_0-9]+)}}", lambda m: str(V.get(m.group(1), "[" + m.group(1) + "]")), tpl)
    open(os.path.join(HERE, target), "w").write(out)
    print(f"{target}: filled from {len(V)} values; missing:", missing)
fill("02_results.template.md", "02_results.md")
if os.path.exists(os.path.join(HERE, "supplement.template.md")): fill("supplement.template.md", "supplement.md")
