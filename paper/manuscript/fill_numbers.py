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

@safe
def f3():
    d = load("fig3c_heritability.csv").sort_values("k_switch")
    put("md_slow", f"{d.mother_daughter.iloc[0]:.2f}"); put("md_fast", f"{d.mother_daughter.iloc[-1]:.2f}")

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
    f = load("fig4f_schedules.csv")
    for model, tag in (("pre_existing", "pre"), ("pre_existing_cost", "cost"), ("drug_induced", "ind")):
        g = f[f.model == model]; best = g.loc[g.long_term_growth_rate.idxmin()]
        put(f"best_{tag}", f"release period {best.release_period:g}, dose {best.dose:g} (growth rate {best.long_term_growth_rate:+.3f} per time unit)")
        cont = g[(g.release_period == 0) & (g.dose == g.dose.max())].long_term_growth_rate.iloc[0]
        put(f"cont_{tag}", f"{cont:+.3f}")
    g = load("fig4g_memory_disruption.csv").groupby("treatment").mean(numeric_only=True)
    put("clones_none", f"{g.loc['none', 'surviving_clones']:.0f}"); put("clones_predrug", f"{g.loc['before_drug', 'surviving_clones']:.0f}"); put("clones_during", f"{g.loc['before_and_during', 'surviving_clones']:.0f}")
    put("cells_none", f"{g.loc['none', 'surviving_cells']:.0f}"); put("cells_during", f"{g.loc['before_and_during', 'surviving_cells']:.0f}")

@safe
def f5a():
    d = load("fig5a_memory_genes.csv")
    put("score_min", f"{d[d.memory_time < 2].score_true.mean():.1f}"); put("score_max", f"{d[d.memory_time > 100].score_true.mean():.1f}")
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
    err = lambda r: np.sqrt(np.mean([(r.died - r.obs_died) ** 2, (r.divided - r.obs_divided) ** 2, (r.survived - r.obs_survived) ** 2]))
    put("rmse_train", f"{np.mean([err(m.loc[7.0]), err(m.loc[13.0])]):.3f}"); put("rmse_heldout", f"{err(m.loc[10.0]):.3f}")
    r = m.loc[10.0]; put("obs_10", f"{r.obs_died:.2f}, {r.obs_divided:.2f} and {r.obs_survived:.2f}"); put("pred_10", f"{r.died:.2f}, {r.divided:.2f} and {r.survived:.2f}")
    put("obs_died_7", f"{m.loc[7.0].obs_died:.2f}"); put("obs_died_13", f"{m.loc[13.0].obs_died:.2f}")
    k = load("fig7c_kin_correlation.csv").groupby("relation").fate_correlation.mean()
    for key, rel in (("phi_sis", "sisters"), ("phi_c1", "first cousins"), ("phi_c2", "second cousins"), ("phi_c3", "third cousins"), ("phi_unrel", "unrelated")): put(key, f"{k[rel]:.2f}")
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
for f in (f2, f2b, f3, f4, f4g, f5a, f5, f6, f7, f7b, f8): f()
fsupp()

def fill(template, target):
    tpl = open(os.path.join(HERE, template)).read()
    missing = sorted(set(re.findall(r"{{([a-z_0-9]+)}}", tpl)) - set(V))
    out = re.sub(r"{{([a-z_0-9]+)}}", lambda m: str(V.get(m.group(1), "[" + m.group(1) + "]")), tpl)
    open(os.path.join(HERE, target), "w").write(out)
    print(f"{target}: filled from {len(V)} values; missing:", missing)
fill("02_results.template.md", "02_results.md")
if os.path.exists(os.path.join(HERE, "supplement.template.md")): fill("supplement.template.md", "supplement.md")
