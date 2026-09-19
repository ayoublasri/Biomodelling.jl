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
    try: fn()
    except Exception as e: print(f"  [{fn.__name__}: {e}]")

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
    put("mem_time", "50")
    d = load("fig4c_decay_vs_dose.csv").sort_values("dose")
    lo = d[d.dose == 0.5].decay_rate.iloc[0]; hi = d[d.dose == 2.0].decay_rate.iloc[0]
    put("decay_lo", f"{lo:.3f}"); put("decay_hi", f"{hi:.3f}"); put("decay_fold", f"{hi / lo:.1f}")
    dd = d[d.dose > 0]
    put("divtime_change", f"{abs(dd.division_time_mean.iloc[-1] - dd.division_time_mean.iloc[0]) / dd.division_time_mean.iloc[0] * 100:.1f}%")
    put("deathtime_range", f"{dd.death_time_mean.max():.1f} to {dd.death_time_mean.min():.1f}")
    c = load("fig4d_fate_concordance.csv"); r = c[(c.model == "memory") & (c.relation == "sisters")].iloc[0]
    put("sis_conc", f"{r.concordance:.2f}"); put("sis_exp", f"{r.expected_independent:.2f}")
    rf = c[(c.model == "fast") & (c.relation == "sisters")].iloc[0]; put("sis_conc_fast", f"{rf.concordance:.2f}")
    e = load("fig4e_clone_diversity.csv"); m = e.groupby("model").mean(numeric_only=True)
    put("clones_before", f"{m.loc['pre_existing', 'effective_clones_before']:.0f}"); put("clones_after_pre", f"{m.loc['pre_existing', 'effective_clones_after']:.0f}")
    put("clones_after_ind", f"{m.loc['drug_induced', 'effective_clones_after']:.0f} of {m.loc['drug_induced', 'effective_clones_before']:.0f}")
    g = load("fig4g_memory_disruption.csv").groupby("treatment").mean(numeric_only=True)
    put("clones_np", f"{g.loc['no_pretreatment', 'surviving_clones']:.0f}"); put("clones_p", f"{g.loc['pretreatment', 'surviving_clones']:.0f}")

@safe
def f5():
    d = load("fig5a_memory_genes.csv")
    put("score_min", f"{d[d.memory_time < 2].score_true.mean():.1f}"); put("score_max", f"{d[d.memory_time > 100].score_true.mean():.1f}")
    y = (d.memory_time > 20).astype(int)
    put("auc_true", f"{roc_auc_score(y, d.score_true):.2f}"); put("auc_seq", f"{roc_auc_score(y, d.score_seq):.2f}")
    m = load("fig5bc_metrics.csv")
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
tpl = open(os.path.join(HERE, "02_results.template.md")).read()
missing = sorted(set(re.findall(r"{{([a-z_0-9]+)}}", tpl)) - set(V))
out = re.sub(r"{{([a-z_0-9]+)}}", lambda m: str(V.get(m.group(1), "[" + m.group(1) + "]")), tpl)
open(os.path.join(HERE, "02_results.md"), "w").write(out)
print("filled", len(V), "placeholders; missing:", missing)
