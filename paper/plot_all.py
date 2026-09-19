"""Generate all figures of the Biomodelling.jl 2.0 preprint from paper/output/*.csv.
Colors: validated colorblind-safe categorical order (blue, orange, aqua, yellow, magenta, green, violet, red)."""
import os, sys, traceback, math
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output"); FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#1c5cab", "#104281"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5, "legend.fontsize": 6.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.labelcolor": INK, "text.color": INK, "legend.frameon": False, "lines.linewidth": 1.4, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
                     "figure.dpi": 150, "savefig.dpi": 300, "pdf.fonttype": 42, "axes.axisbelow": True})

def load(name):
    p = os.path.join(OUT, name)
    if not os.path.exists(p): raise FileNotFoundError(name)
    return pd.read_csv(p)
def kv(name):
    d = load(name); return dict(zip(d.key, d.value))
def label(ax, letter, title=None):
    ax.text(-0.18, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    if title: ax.set_title(title, loc="left", fontsize=7.5, color=INK2)
def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight"); fig.savefig(os.path.join(FIG, name + ".png"), bbox_inches="tight"); plt.close(fig)
    print("saved", name)
def panel(fn):
    def wrapped(*a, **k):
        try: fn(*a, **k)
        except Exception as e: print(f"  [skipped {fn.__name__}: {e}]")
    return wrapped

# ------------------------------------------------------------------ Figure 2
def fig2():
    fig = plt.figure(figsize=(7.2, 4.6)); gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.55, wspace=0.45)
    @panel
    def a(ax):
        d = load("fig2a_birthdeath.csv")
        ax.bar(d.n, d.empirical, color=C[0], width=0.8, label="simulation")
        ax.plot(d.n, d.poisson, color=INK, lw=1.2, marker="o", ms=2.5, label="Poisson(12)")
        ax.set_xlabel("molecules"); ax.set_ylabel("probability"); ax.legend(); label(ax, "a", "birth–death")
    @panel
    def b(ax):
        d = load("fig2b_telegraph.csv")
        for i, (reg, g) in enumerate(d.groupby("regime", sort=False)):
            ax.bar(g.n, g.empirical, color=C[i], width=0.8, alpha=0.75, label=f"simulation ({reg})")
            ax.plot(g.n, g.beta_poisson, color=INK if i == 0 else INK2, lw=1.0, ls="-" if i == 0 else "--", label=f"Beta-Poisson ({reg})")
        ax.set_xlim(-0.5, 60); ax.set_xlabel("mRNA molecules"); ax.set_ylabel("probability"); ax.legend(); label(ax, "b", "telegraph model")
    @panel
    def c(ax):
        d = load("fig2c_bursty.csv")
        ax.bar(d.n, d.empirical, color=C[2], width=0.8, label="simulation")
        ax.plot(d.n, d.negbin, color=INK, lw=1.2, label="negative binomial")
        ax.set_xlabel("protein molecules"); ax.set_ylabel("probability"); ax.legend(); label(ax, "c", "bursty protein")
    @panel
    def d(ax):
        d = load("fig2d_kernels.csv")
        y = np.arange(len(d)); h = 0.38
        ax.barh(y - h/2, d.ks_telegraph, height=h, color=C[0], label="telegraph")
        ax.barh(y + h/2, d.ks_birthdeath, height=h, color=C[1], label="birth–death")
        crit = 1.63 / math.sqrt(20000)
        ax.axvline(crit, color=INK2, ls=":", lw=0.9); ax.text(crit, len(d) - 0.4, " KS 99% critical", color=INK2, fontsize=6, va="top")
        ax.set_yticks(y); ax.set_yticklabels(d.kernel); ax.invert_yaxis(); ax.set_xlabel("KS distance to exact law (n = 20 000)"); ax.legend(loc="lower right"); label(ax, "d", "kernels")
        for i, r in d.iterrows(): ax.text(max(r.ks_telegraph, r.ks_birthdeath) + 0.0005, i, f"{r.sec_per_1e4_telegraph:.2f} s / 10⁴ cells", va="center", fontsize=6, color=INK2)
    @panel
    def e(ax):
        d = load("fig2e_crosscheck.csv"); s = kv("fig2e_crosscheck_stats.csv")
        ax.plot(d.n, d.exact, color=INK, lw=1.0, label="exact")
        ax.plot(d.n, d.biomodelling, color=C[0], lw=0, marker="o", ms=3, label="Biomodelling.jl")
        ax.plot(d.n, d.jumpprocesses, color=C[1], lw=0, marker="s", ms=2.5, mfc="none", label="JumpProcesses.jl")
        ax.set_xlim(-0.5, 50); ax.set_xlabel("mRNA molecules"); ax.set_ylabel("probability"); ax.legend()
        ax.text(0.98, 0.55, f"two-sample KS = {s['ks_two_sample']:.4f}\n(99% critical {s['ks_two_sample_crit99']:.4f})\n{s['time_biomodelling_serial_s']:.1f} s vs {s['time_jumpprocesses_serial_s']:.1f} s", transform=ax.transAxes, ha="right", va="top", fontsize=6, color=INK2)
        label(ax, "e", "independent implementation")
    @panel
    def f(ax):
        d = load("fig2f_runtime.csv")
        for i, (k, g) in enumerate(d.groupby("kernel", sort=False)):
            g10 = g[g.genes == 10].sort_values("cells")
            ax.plot(g10.cells, g10.seconds, marker="o", ms=3, color=C[i], label=f"{k}, 10 genes")
            g1000 = g[g.cells == 1000].sort_values("genes")
            ax.plot(g1000.genes * 100, g1000.seconds, marker="s", ms=3, ls="--", color=C[i], label=f"{k}, 1000 cells (×100 genes)")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("cells  (or 100 × genes)"); ax.set_ylabel("seconds for 200 steps"); ax.legend(fontsize=5.5); label(ax, "f", "runtime")
    for fn, pos in zip((a, b, c, d, e, f), range(6)):
        fn(fig.add_subplot(gs[pos // 3, pos % 3]))
    save(fig, "fig2_engine")

# ------------------------------------------------------------------ Figure 3
def fig3():
    fig = plt.figure(figsize=(7.2, 5.0)); gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.6, wspace=0.45)
    @panel
    def a(cell):
        d = load("fig3a_trace.csv")
        sub = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=cell, hspace=0.15)
        axes = [fig.add_subplot(sub[i]) for i in range(3)]
        on = d.G_on.values > 0
        for ax in axes:
            ax.fill_between(d.t, 0, 1, where=on, transform=ax.get_xaxis_transform(), color=C[3], alpha=0.18, lw=0)
        axes[0].plot(d.t, d.V, color=C[0]); axes[0].set_ylabel("V")
        axes[1].plot(d.t, d.mRNA, color=C[1]); axes[1].set_ylabel("mRNA")
        axes[2].plot(d.t, d.protein, color=C[2]); axes[2].set_ylabel("protein"); axes[2].set_xlabel("time")
        for ax in axes[:2]: ax.set_xticklabels([])
        label(axes[0], "a", "single lineage (shading: promoter on)")
    @panel
    def b(ax):
        d0 = load("fig3b_scaling_no_replication.csv"); d1 = load("fig3b_scaling_replication.csv")
        ax.scatter(d0.volume, d0.mRNA, s=3, color=C[0], alpha=0.35, lw=0, label="no replication")
        conc = (d0.mRNA / d0.volume).mean(); xs = np.linspace(d0.volume.min(), d0.volume.max(), 10)
        ax.plot(xs, conc * xs, color=INK, lw=1.0, label=f"mean concentration × V")
        bins = np.linspace(1, 2, 9); idx = np.digitize(d1.volume, bins)
        m = [d1.mRNA[idx == i].mean() for i in range(1, len(bins))]; ctr = 0.5 * (bins[1:] + bins[:-1])
        ax.plot(ctr, m, color=C[1], marker="o", ms=3, label="binned mean, with replication")
        ax.set_xlabel("cell volume"); ax.set_ylabel("mRNA molecules"); ax.legend(fontsize=5.5); label(ax, "b", "size scaling")
    @panel
    def c(ax):
        d = load("fig3c_heritability.csv")
        for i, (col, lab) in enumerate((("mother_daughter", "mother–daughter"), ("sisters", "sisters"), ("cousins", "cousins"))):
            ax.plot(d.k_switch, d[col], marker="o", ms=3, color=C[i], label=lab)
        ax.plot(d.k_switch, d.mother_daughter_mRNA, ls="--", color=C[0], lw=1.0, label="mother–daughter (mRNA)")
        ax.axvline(math.log(2) / 20, color=INK2, ls=":", lw=0.8); ax.text(math.log(2) / 20, 0.95, " 1/cycle time", fontsize=6, color=INK2, va="top")
        ax.set_xscale("log"); ax.set_xlabel("promoter switching rate k_on = k_off"); ax.set_ylabel("correlation (protein at division)"); ax.set_ylim(-0.1, 1); ax.legend(fontsize=5.5, loc="lower left"); label(ax, "c", "heritability")
    @panel
    def d(ax):
        d = load("fig3d_memory_timescale.csv")
        for i, (k, g) in enumerate(d.groupby("k_switch")):
            g = g.sort_values("cycle_time_nominal")
            ax.plot(g.cycle_time_measured, g.tau_time, marker="o", ms=3.5, color=C[i], label=f"k = {k:g}  (1/2k = {1/(2*k):.0f})")
            ax.axhline(1 / (2 * k), color=C[i], ls=":", lw=0.8)
        ax.set_ylim(0, None); ax.set_xlabel("cell-cycle time"); ax.set_ylabel("memory timescale (time units)"); ax.legend(fontsize=5.5); label(ax, "d", "memory timescale vs cell-cycle time")
    @panel
    def e(ax):
        d = load("fig3e_noise.csv"); d = d[d.partition_sigma == d.partition_sigma.max()]
        cvs = sorted(d.cycle_cv.unique()); x = np.arange(len(cvs)); w = 0.2
        for i, (part, hatch) in enumerate((("binomial", ""), ("betabinomial", "///"))):
            g = d[d.partitioning == part].set_index("cycle_cv").loc[cvs]
            ax.bar(x + (2*i - 1.5) * w, g.cv2_population, width=w, color=C[0], hatch=hatch, edgecolor="white", label=f"population, {part}")
            ax.bar(x + (2*i - 0.5) * w, g.cv2_lineage, width=w, color=C[1], hatch=hatch, edgecolor="white", label=f"single lineage, {part}")
        ax.set_xticks(x); ax.set_xticklabels([f"{c:g}" for c in cvs]); ax.set_xlabel("cell-cycle time CV (partition σ = 0.1)")
        ax.set_ylabel("CV² of protein concentration"); ax.legend(fontsize=5); label(ax, "e", "population vs lineage noise")
    @panel
    def f(ax):
        d = load("fig3f_copynumber.csv")
        x = np.arange(3); w = 0.26
        series = (("k_on", d.k_on.values, 0.5, "k_on (truth 0.5)"), ("k_off", d.k_off.values, 1.5, "k_off (truth 1.5)"),
                  ("k_tx", d.k_tx.values / d.mean_volume.values / 20, 2.0, "k_tx / (V̄ · 20)  (truth 2 per allele)"))
        for i, (col, vals, tv, lab) in enumerate(series):
            ax.bar(x + (i - 1) * w, vals, width=w, color=C[i], label=lab)
            ax.hlines(tv, x[0] - 0.45, x[-1] + 0.45, colors=C[i], linestyles=":", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(["1 copy", "2 copies", "pooled"]); ax.set_ylabel("fitted telegraph parameters"); ax.legend(fontsize=5); label(ax, "f", "copy number and naive fits")
    a(gs[0, 0]); b(fig.add_subplot(gs[0, 1])); c(fig.add_subplot(gs[0, 2])); d(fig.add_subplot(gs[1, 0])); e(fig.add_subplot(gs[1, 1])); f(fig.add_subplot(gs[1, 2]))
    save(fig, "fig3_growth")

# ------------------------------------------------------------------ Figure 4
def fig4():
    fig = plt.figure(figsize=(7.2, 7.0)); gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.6, wspace=0.5)
    @panel
    def a(ax):
        d = load("fig4a_trajectories.csv")
        for i, (s, g) in enumerate(d.groupby("schedule", sort=False)):
            ax.plot(g.t, g.popsize, color=C[i], label=s.replace("_", " "))
        g = d[d.schedule == d.schedule.iloc[0]]
        ax.fill_between(g.t, 1, 1e9, where=g.dose > 0, color=INK2, alpha=0.08, lw=0)
        ax.set_yscale("log"); ax.set_ylim(d.popsize.min() * 0.5, d.popsize.max() * 2); ax.set_xlabel("time"); ax.set_ylabel("population size"); ax.legend(fontsize=5.5); label(ax, "a", "schedules")
    @panel
    def b(ax):
        d = load("fig4b_killcurves.csv")
        for i, (dose, g) in enumerate(d.groupby("dose")):
            ax.plot(g.t_since_drug, g.surviving_fraction, color=C[i], label=f"dose {dose:g}")
        ax.set_yscale("log"); ax.set_xlabel("time since drug"); ax.set_ylabel("N(t) / N(0)"); ax.legend(fontsize=5.5); label(ax, "b", "kill curves")
    @panel
    def c(ax):
        d = load("fig4c_decay_vs_dose.csv")
        ax.plot(d.dose, d.decay_rate, marker="o", ms=4, color=C[0]); ax.set_xlabel("dose"); ax.set_ylabel("population decay rate", color=C[0])
        ax2 = ax.inset_axes([0.55, 0.55, 0.42, 0.4])
        ax2.plot(d.dose, d.division_time_mean, marker="s", ms=3, color=C[1]); ax2.set_title("division time", fontsize=6, color=INK2); ax2.tick_params(labelsize=5); ax2.set_ylim(0, d.division_time_mean.max() * 1.3)
        label(ax, "c", "decay vs single-cell timing")
    @panel
    def d(ax):
        d = load("fig4d_fate_concordance.csv")
        x = np.arange(len(d)); w = 0.38
        ax.bar(x - w/2, d.concordance, width=w, color=C[0], label="observed")
        ax.bar(x + w/2, d.expected_independent, width=w, color=INK2, alpha=0.5, label="independent fates")
        ax.set_xticks(x); ax.set_xticklabels([f"{r.model}\n{r.relation}" for _, r in d.iterrows()], fontsize=6); ax.set_ylabel("fate concordance"); ax.set_ylim(0, 1); ax.legend(fontsize=5.5); label(ax, "d", "related cells share fates")
    @panel
    def e(ax):
        d = load("fig4e_clone_diversity.csv"); m = d.groupby("model").agg(["mean", "std"])
        x = np.arange(len(m)); w = 0.38
        ax.bar(x - w/2, m["effective_clones_before"]["mean"], yerr=m["effective_clones_before"]["std"], width=w, color=C[0], label="before drug", capsize=2)
        ax.bar(x + w/2, m["effective_clones_after"]["mean"], yerr=m["effective_clones_after"]["std"], width=w, color=C[1], label="after drug", capsize=2)
        ax.set_xticks(x); ax.set_xticklabels([i.replace("_", " ") for i in m.index]); ax.set_ylabel("effective number of clones"); ax.legend(fontsize=5.5); label(ax, "e", "barcode diversity")
    @panel
    def f(cell):
        d = load("fig4f_schedules.csv")
        sub = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=cell, wspace=0.15)
        vmax = np.abs(d.long_term_growth_rate).max()
        for i, model in enumerate(("pre_existing", "drug_induced")):
            ax = fig.add_subplot(sub[i]); g = d[d.model == model]
            piv = g.pivot(index="release_period", columns="dose", values="long_term_growth_rate")
            im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto", origin="lower")
            ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels([f"{c:g}" for c in piv.columns]); ax.set_yticks(range(len(piv.index))); ax.set_yticklabels([f"{r:g}" for r in piv.index])
            ax.set_xlabel("dose"); ax.set_title(model.replace("_", " "), fontsize=7, color=INK2); ax.grid(False)
            if i == 0: ax.set_ylabel("release period"); label(ax, "f", None)
            else: ax.set_yticklabels([])
            for (r, c), v in np.ndenumerate(piv.values): ax.text(c, r, f"{v:+.3f}", ha="center", va="center", fontsize=4.5, color=INK)
        fig.colorbar(im, ax=ax, fraction=0.08, pad=0.04, label="long-term growth rate")
    @panel
    def g(ax):
        d = load("fig4g_memory_disruption.csv"); m = d.groupby("treatment", sort=False).agg(["mean", "std"])
        x = np.arange(len(m))
        ax.bar(x, m["surviving_clones"]["mean"], yerr=m["surviving_clones"]["std"], color=[C[0], C[1]], width=0.6, capsize=2)
        for i, t in enumerate(m.index): ax.text(i, m["surviving_clones"]["mean"].iloc[i] + m["surviving_clones"]["std"].iloc[i] + 1, f"{m['surviving_clones']['mean'].iloc[i]:.0f}", ha="center", fontsize=6)
        ax.set_xticks(x); ax.set_xticklabels([t.replace("_", " ") for t in m.index]); ax.set_ylabel("surviving clones (resistant colonies)"); label(ax, "g", "memory disruption before drug")
    @panel
    def h(ax):
        d = load("fig4h_mgmt.csv")
        for i, (m, g) in enumerate(d.groupby("model", sort=False)):
            ax.plot(g.t, g.high_fraction, color=C[i], label=m.replace("_", " "))
        g = d[d.model == d.model.iloc[0]]; ax.fill_between(g.t, 0, 1, where=g.dose > 0, color=INK2, alpha=0.08, lw=0)
        ax.set_xlabel("time"); ax.set_ylabel("fraction of high-expressing cells"); ax.set_ylim(0, 1); ax.legend(fontsize=5.5); label(ax, "h", "phenotypic selection persists with slow switching")
    a(fig.add_subplot(gs[0, 0])); b(fig.add_subplot(gs[0, 1])); c(fig.add_subplot(gs[0, 2])); d(fig.add_subplot(gs[1, 0])); e(fig.add_subplot(gs[1, 1])); f(gs[1, 2]); g(fig.add_subplot(gs[2, 0])); h(fig.add_subplot(gs[2, 1]))
    save(fig, "fig4_persisters")

# ------------------------------------------------------------------ Figure 5
def fig5():
    fig = plt.figure(figsize=(7.2, 4.8)); gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.4)
    @panel
    def a(ax):
        d = load("fig5a_memory_genes.csv")
        ax.plot(d.memory_time, d.score_true, marker="o", ms=3, lw=0, color=C[0], label="true concentrations")
        ax.plot(d.memory_time, d.score_seq, marker="s", ms=3, lw=0, mfc="none", color=C[1], label="sequenced counts")
        ax.axvline(20, color=INK2, ls=":", lw=0.8); ax.text(20, ax.get_ylim()[1] if ax.get_ylim()[1] > 1 else 2, " cycle time", fontsize=6, color=INK2)
        ax.axhline(1, color=INK2, lw=0.6); ax.set_xscale("log"); ax.set_xlabel("memory time 1/(k_on+k_off)"); ax.set_ylabel("clonal variance score"); ax.legend(fontsize=5.5); label(ax, "a", "memory genes are recoverable")
    @panel
    def b(ax):
        d = load("fig5bc_metrics.csv"); d = d[~d.dataset.str.startswith("imputed")]
        order = ["fixed_volume", "population_counts", "population_concentration", "sequenced_counts", "sequenced_normalized"]
        methods = [m for m in ["pearson", "spearman", "genie3"] if m in set(d.method)]
        x = np.arange(len(order)); w = 0.8 / len(methods)
        for i, m in enumerate(methods):
            vals = [d[(d.method == m) & (d.dataset == o)].aupr.mean() if ((d.method == m) & (d.dataset == o)).any() else np.nan for o in order]
            ax.bar(x + (i - (len(methods) - 1) / 2) * w, vals, width=w, color=C[i], label=m)
        ax.axhline(d.random_aupr.mean(), color=INK2, ls=":", lw=0.8); ax.text(len(order) - 0.5, d.random_aupr.mean(), "random", fontsize=6, color=INK2, ha="right", va="bottom")
        ax.set_xticks(x); ax.set_xticklabels([o.replace("_", "\n") for o in order], fontsize=5.5); ax.set_ylabel("AUPR"); ax.legend(fontsize=5.5); label(ax, "b", "GRN inference under growth and division")
    @panel
    def c(ax):
        d = load("fig5bc_metrics.csv")
        order = ["sequenced_counts", "imputed_knn_smoothing", "imputed_magic"]; order = [o for o in order if o in set(d.dataset)]
        methods = [m for m in ["pearson", "genie3"] if m in set(d.method)]
        x = np.arange(len(order)); w = 0.8 / len(methods)
        for i, m in enumerate(methods):
            vals = [d[(d.method == m) & (d.dataset == o)].aupr.mean() for o in order]
            ax.bar(x + (i - (len(methods) - 1) / 2) * w, vals, width=w, color=C[i], label=m)
        ax.set_xticks(x); ax.set_xticklabels([o.replace("imputed_", "").replace("_", "\n") for o in order], fontsize=5.5); ax.set_ylabel("AUPR"); ax.legend(fontsize=5.5); label(ax, "c", "imputation before network inference")
    @panel
    def d(ax):
        d = load("fig5d_summary.csv")
        x = np.arange(len(d)); w = 0.38
        ax.bar(x - w/2, d.r2_zero_baseline, width=w, color=INK2, alpha=0.5, label="no-change baseline")
        ax.bar(x + w/2, d.r2_correlation_baseline, width=w, color=C[0], label="correlation propagation")
        ax.axhline(0, color=INK, lw=0.6)
        ax.set_xticks(x); ax.set_xticklabels([f"KD gene {int(g)}" for g in d.knocked_gene], fontsize=6); ax.set_ylabel("R² of predicted log₂ fold changes"); ax.legend(fontsize=5.5); label(ax, "d", "perturbation ground truth")
    a(fig.add_subplot(gs[0, 0])); b(fig.add_subplot(gs[0, 1])); c(fig.add_subplot(gs[1, 0])); d(fig.add_subplot(gs[1, 1]))
    save(fig, "fig5_benchmarks")

# ------------------------------------------------------------------ Figure 6
def fig6():
    fig = plt.figure(figsize=(7.2, 4.8)); gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.6, wspace=0.45)
    def posterior_panels(cell, particles_file, summary_file, letter, title, ref_prefix, ref_label):
        p = load(particles_file); s = kv(summary_file)
        sub = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=cell, wspace=0.5)
        for i, k in enumerate(("k_on", "k_off", "k_tx")):
            ax = fig.add_subplot(sub[i])
            ax.hist(np.log10(p[k]), bins=25, weights=p.weight, color=C[0], alpha=0.8)
            ax.axvline(np.log10(s[f"true_{k}"]), color=INK, lw=1.2, label="truth")
            ax.axvline(np.log10(s[f"{ref_prefix}_{k}"]), color=C[7], lw=1.2, ls="--", label=ref_label)
            ax.set_xlabel(f"log₁₀ {k}"); ax.set_yticks([])
            if i == 0: label(ax, letter, title); ax.set_ylabel("posterior")
            if i == 2: ax.legend(fontsize=5.5)
    @panel
    def a(cell): posterior_panels(cell, "fig6a_particles.csv", "fig6a_summary.csv", "a", "non-dividing cells: ABC vs exact MLE", "mle", "exact MLE")
    @panel
    def b(cell): posterior_panels(cell, "fig6b_particles.csv", "fig6b_summary.csv", "b", "dividing population: division-aware ABC vs naive fit", "naive", "naive telegraph fit")
    @panel
    def c(ax):
        p = load("fig6c_particles.csv"); s = kv("fig6c_summary.csv")
        ax.scatter(p.h_max, p.K, s=8 + 400 * p.weight, color=C[0], alpha=0.5, lw=0, label="posterior particles")
        ax.plot([s["true_h_max"]], [s["true_K"]], marker="*", ms=11, color=C[7], lw=0, label="truth")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("h_max"); ax.set_ylabel("K (protection threshold)"); ax.legend(fontsize=5.5); label(ax, "c", "drug parameters from kill curve + sister fates")
    @panel
    def d(ax):
        d = load("fig6d_ppc.csv")
        for c in [c for c in d.columns if c.startswith("replicate")]:
            ax.plot(d.t, d[c], color=C[0], alpha=0.5, lw=1.0)
        ax.plot(d.t, d.observed, color=INK, lw=1.6, label="observed")
        ax.plot([], [], color=C[0], alpha=0.5, label="posterior predictive")
        ax.set_yscale("log"); ax.set_xlabel("time"); ax.set_ylabel("N(t)/N(drug start)"); ax.legend(fontsize=5.5); label(ax, "d", "posterior predictive check")
    a(gs[0, :]); b(gs[1, :2]); c(fig.add_subplot(gs[1, 2]))
    save(fig, "fig6_inference")
    fig2_ = plt.figure(figsize=(3.4, 2.4)); d(fig2_.add_subplot(111)); save(fig2_, "fig6d_ppc")

if __name__ == "__main__":
    which = sys.argv[1:] or ["2", "3", "4", "5", "6"]
    for w in which:
        try: globals()[f"fig{w}"]()
        except Exception: traceback.print_exc()
