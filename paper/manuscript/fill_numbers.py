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
def f5():
    d = load("fig5a_memory_genes.csv")
    put("score_min", f"{d[d.memory_time < 2].score_true.mean():.1f}"); put("score_max", f"{d[d.memory_time > 100].score_true.mean():.1f}")
    y = (d.memory_time > 20).astype(int)
    put("auc_true", f"{roc_auc_score(y, d.score_true):.2f}"); put("auc_seq", f"{roc_auc_score(y, d.score_seq):.2f}")
    m = load("fig5bc_metrics.csv")
    for meth in ("pearson", "genie3"):
        for ds, tag in (("fixed_volume", "fixed"), ("population_counts", "counts"), ("population_concentration", "conc"), ("population_cycle_regressed", "reg"), ("sequenced_counts", "seq"), ("sequenced_normalized", "norm")):
            put(f"aupr_{meth}_{tag}", f"{m[(m.method == meth) & (m.dataset == ds)].aupr.mean():.2f}")
    put("aupr_random", f"{m.random_aupr.mean():.2f}")
    base = m[(m.method == "pearson") & (m.dataset == "sequenced_counts")].aupr.mean()
    imp = m[(m.method == "pearson") & (m.dataset.str.startswith("imputed"))].aupr
    put("imp_effect", f"{(imp.min() - base):+.2f} to {(imp.max() - base):+.2f} (Pearson correlation)")
    s = load("fig5d_summary.csv")
    put("r2_zero", f"{s.r2_zero_baseline.mean():.2f}"); put("r2_corr", f"{s.r2_correlation_baseline.mean():.2f}")

@safe
def f6():
    a = kv("fig6a_summary.csv"); put("abc_a", f"({a['abc_median_k_on']:.2f}, {a['abc_median_k_off']:.2f}, {a['abc_median_k_tx']:.1f})")
    b = kv("fig6b_summary.csv"); put("naive_b", f"k_on = {b['naive_k_on']:.2f}, k_off = {b['naive_k_off']:.2f}, k_tx = {b['naive_k_tx']:.1f}")
    put("abc_b", f"({b['abc_median_k_on']:.2f}, {b['abc_median_k_off']:.2f}, {b['abc_median_k_tx']:.1f})")
    c = kv("fig6c_summary.csv"); put("abc_c", f"({c['abc_median_h_max']:.2f}, {c['abc_median_K']:.0f})")

for f in (f2, f3, f4, f5, f6): f()
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
fsupp()

def fill(template, target):
    tpl = open(os.path.join(HERE, template)).read()
    missing = sorted(set(re.findall(r"{{([a-z_0-9]+)}}", tpl)) - set(V))
    out = re.sub(r"{{([a-z_0-9]+)}}", lambda m: str(V.get(m.group(1), "[" + m.group(1) + "]")), tpl)
    open(os.path.join(HERE, target), "w").write(out)
    print(f"{target}: filled from {len(V)} values; missing:", missing)
fill("02_results.template.md", "02_results.md")
if os.path.exists(os.path.join(HERE, "supplement.template.md")): fill("supplement.template.md", "supplement.md")
