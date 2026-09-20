"""Supplementary figures: S1 death/division-time distributions per dose, S2 ABC tolerance schedules, S3 runtime table."""
import os, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "output"); FIG = os.path.join(HERE, "figures")
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
try:                                    # the typeface and spacing of the main figures
    from matplotlib import font_manager
    for w in ("400", "600", "800"):
        fp = os.path.join(FIG, "fonts", f"Inter-{w}.ttf")
        if os.path.exists(fp): font_manager.fontManager.addfont(fp)
    FAMILY = "Inter" if any(f.name == "Inter" for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
except Exception:
    FAMILY = "DejaVu Sans"
plt.rcParams.update({"font.family": FAMILY, "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
                     "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.7,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.labelcolor": INK, "text.color": INK, "axes.titlepad": 9, "axes.labelpad": 3.5,
                     "legend.frameon": False, "legend.handlelength": 1.5, "legend.labelspacing": 0.35,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.45, "axes.axisbelow": True,
                     "savefig.dpi": 400, "pdf.fonttype": 42, "figure.facecolor": "white", "savefig.facecolor": "white"})
def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight"); fig.savefig(os.path.join(FIG, name + ".png"), bbox_inches="tight"); plt.close(fig); print("saved", name)
# S1: single-cell timing distributions
try:
    d = pd.read_csv(os.path.join(OUT, "fig4c_times.csv"))
    dcol = {dose: C[i % 8] for i, dose in enumerate(sorted(d.dose.unique()))}  # one colour per dose across panels (matches Fig. 4b)
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
    for dose, g in d[d.event == "division"].groupby("dose"):
        axes[0].hist(g.time, bins=np.linspace(10, 35, 40), histtype="step", color=dcol[dose], label=f"dose {dose:g}", density=True)
    axes[0].set_xlabel("division time (cells born under drug)"); axes[0].set_ylabel("density"); axes[0].legend(fontsize=6); axes[0].set_title("a   division times", loc="left", fontweight="semibold", color="#2a2a28")
    death = d[d.event == "death"]
    for dose, g in death.groupby("dose"):
        axes[1].hist(g.time, bins=np.linspace(0, 60, 40), histtype="step", color=dcol[dose], label=f"dose {dose:g}", density=True)
    axes[1].set_xlim(0, float(death.time.quantile(0.99)) * 1.1); axes[1].set_xlabel("time from drug (or birth) to death"); axes[1].legend(fontsize=6); axes[1].set_title("b   death times", loc="left", fontweight="semibold", color="#2a2a28")
    save(fig, "figS1_timing")
except Exception as e: print("S1 skipped:", e)
# S2: ABC schedules
try:
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
    for i, (tag, lab) in enumerate((("fig6a", "telegraph, non-dividing"), ("fig6b", "telegraph, dividing population"), ("fig6c", "drug parameters"))):
        s = pd.read_csv(os.path.join(OUT, f"{tag}_schedule.csv"))
        axes[0].plot(s.generation[1:], s.epsilon[1:], marker="o", ms=3, color=C[i], label=lab)
        axes[1].plot(s.generation[1:], s.acceptance[1:], marker="o", ms=3, color=C[i], label=lab)
    axes[0].set_yscale("log"); axes[0].set_xlabel("generation"); axes[0].set_ylabel("tolerance ε"); axes[0].legend(fontsize=6); axes[0].set_title("a  tolerance schedule", loc="left")
    axes[1].set_yscale("log"); axes[1].set_xlabel("generation"); axes[1].set_ylabel("acceptance rate"); axes[1].set_title("b  acceptance rate", loc="left")
    save(fig, "figS2_abc_schedules")
except Exception as e: print("S2 skipped:", e)
# S3: runtime table (markdown)
try:
    r = pd.read_csv(os.path.join(OUT, "fig2f_runtime.csv"))
    with open(os.path.join(HERE, "manuscript", "tableS_runtime.md"), "w") as f:
        f.write("| genes | cells | kernel | threads | seconds (200 steps) | cell-steps / s |\n|---|---|---|---|---|---|\n")
        for _, x in r.iterrows(): f.write(f"| {int(x.genes)} | {int(x.cells)} | {x.kernel} | {int(x.threads)} | {x.seconds:.2f} | {x.cell_steps_per_second:,.0f} |\n")
    print("runtime table written")
except Exception as e: print("S3 skipped:", e)

# S3: identifiability and numerical robustness of the calibration
try:
    import numpy as _np
    prof = pd.read_csv(os.path.join(OUT, "fig7g_profile.csv"))
    slic = pd.read_csv(os.path.join(OUT, "fig7h_memory_slice.csv"))
    step = pd.read_csv(os.path.join(OUT, "fig7i_stepsize.csv"))
    cal = pd.read_csv(os.path.join(OUT, "fig7_calibration.csv")).set_index("key").value
    obs = pd.read_csv(os.path.join(HERE, "data", "iyer2025_u2os_fates.csv"))
    ses = []
    for _, r in obs.iterrows():
        n = r["cells_at_drug"]
        for k in ("died", "divided", "survived_without_dividing"):
            q = r[k] / n; ses.append(_np.sqrt(q * (1 - q) / n))
    se = float(_np.mean(ses))

    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 0.72], wspace=0.72)

    ax = fig.add_subplot(gs[0])
    best = prof.rmse.min()
    pal6 = [C[0], C[1], C[2], C[3], C[4], C[7]]  # avoid two near-greens (IC50 -> C[7])
    plabels = {"k_on": "k_on", "k_off": "k_off", "h_max": "h_max", "EC50": "EC50", "m_h": "Hill m", "IC50": "IC50 growth"}
    for i, (par, g) in enumerate(prof.groupby("parameter", sort=False)):
        g = g.sort_values("value")
        ax.plot(g.value / float(cal[par]), g.rmse, marker="o", ms=2.2, lw=1.1, color=pal6[i % len(pal6)], label=plabels.get(par, par))
    ax.axhline(best + se, color="#52514e", ls=":", lw=0.9, label="within one s.e. of the data")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("parameter / fitted value"); ax.set_ylabel("RMSE of fate fractions")
    ax.legend(fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.30), columnspacing=1.2, handlelength=1.2)
    ax.set_title("a   one parameter at a time", loc="left", fontweight="semibold", color="#2a2a28")

    ax = fig.add_subplot(gs[1])
    piv = slic.pivot(index="p_on", columns="memory_generations", values="rmse")
    im = ax.imshow(piv.values, origin="lower", aspect="auto", cmap="viridis_r",
                   extent=[-0.5, piv.shape[1] - 0.5, -0.5, piv.shape[0] - 0.5])
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if piv.values[i, j] <= slic.rmse.min() + se:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#e34948", lw=1.0))
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels([f"{c:g}" for c in piv.columns], fontsize=6)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels([f"{r:g}" for r in piv.index], fontsize=6)
    ax.set_xlabel("memory (generations)"); ax.set_ylabel("fraction resistant")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02); cb.set_label("RMSE", fontsize=6.5, labelpad=1); cb.ax.tick_params(labelsize=6)
    ax.set_title("b   memory and resistant fraction", loc="left", fontweight="semibold", color="#2a2a28")

    ax = fig.add_subplot(gs[2])
    g = step.groupby("dt_h")[["died", "divided", "survived"]].agg(["mean", "std"])
    x = _np.arange(3); w = 0.36
    for k, dt in enumerate(sorted(step.dt_h.unique())):
        ax.bar(x + (k - 0.5) * w, [g.loc[dt, (c, "mean")] for c in ("died", "divided", "survived")],
               yerr=[g.loc[dt, (c, "std")] for c in ("died", "divided", "survived")],
               width=w, color=C[k], capsize=1.5, label=f"dt = {dt:g} h")
    ax.set_xticks(x); ax.set_xticklabels(["died", "divided", "survived"], fontsize=6, rotation=20, ha="right", rotation_mode="anchor")
    ax.set_ylabel("fraction of cells at 13 µM")
    ax.legend(fontsize=6, labelspacing=0.25, loc="upper right")
    ax.set_title("c   integration step", loc="left", fontweight="semibold", color="#2a2a28")
    save(fig, "figS3_identifiability")
except Exception as e:
    print("S3 skipped:", e)

# S4: does cell-cycle-dependent killing explain the cisplatin fates better than heritable expression?
try:
    import numpy as _np
    prof = pd.read_csv(os.path.join(OUT, "fig7j_cycle_profile.csv"))
    fat = pd.read_csv(os.path.join(OUT, "fig7j_cycle_fates.csv"))
    fit = pd.read_csv(os.path.join(OUT, "fig7j_cycle_fit.csv")).set_index("key").value
    obs = pd.read_csv(os.path.join(HERE, "data", "iyer2025_u2os_fates.csv"))
    ses, ov = [], {}
    for _, r in obs.iterrows():
        n = r["cells_at_drug"]
        for k, lab in (("died", "died"), ("divided", "divided"), ("survived_without_dividing", "survived")):
            q = r[k] / n; ses.append(_np.sqrt(q * (1 - q) / n)); ov[(r["cisplatin_uM"], lab)] = q
    se = float(_np.mean(ses))
    held = 10.0

    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.7), gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.42})

    ax = axes[0]
    ax.plot(prof.cycle_baseline, prof.rmse, marker="o", ms=3, lw=1.2, color=C[0])
    b = prof.rmse.min()
    ax.axhline(b + se, color="#52514e", ls=":", lw=0.9, label="within one s.e. of the data")
    ax.axvline(float(fit["cycle_baseline"]), color=C[3], ls="--", lw=1.0, label="fitted value")
    ax.set_xlabel("cycle-independent fraction of the hazard"); ax.set_ylabel("RMSE of fate fractions")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.30))
    ax.set_title("a   is the cycle dependence determined?", loc="left", fontweight="semibold", color="#2a2a28")

    ax = axes[1]
    for k, (mod, col) in enumerate((("cycle-independent", C[0]), ("cycle-dependent", C[1]))):
        g = fat[fat.model == mod].groupby(["cisplatin_uM"])[["died", "divided", "survived"]].mean()
        for dd, row in g.iterrows():
            for lab in ("died", "divided", "survived"):
                o = ov.get((dd, lab))
                if o is None: continue
                ax.plot([o], [row[lab]], marker=("o" if dd != held else "s"), ms=4.5, lw=0,
                        mfc=(col if dd != held else "none"), mec=col, mew=1.0)
        ax.plot([], [], marker="o", ms=4.5, lw=0, color=col, label=mod)
    lim = [0, 0.75]
    ax.plot(lim, lim, color="#52514e", lw=0.8, ls="-")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("observed fraction"); ax.set_ylabel("simulated fraction")
    ax.plot([], [], marker="s", ms=4.5, lw=0, mfc="none", mec="#52514e", label="held-out 10 µM")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2, columnspacing=1.0)
    ax.set_title("b   fates under the two mechanisms", loc="left", fontweight="semibold", color="#2a2a28")
    save(fig, "figS4_cycle")
except Exception as e:
    print("S4 skipped:", e)

# S5: do the persister and melanoma conclusions survive cell-cycle-dependent killing?
try:
    import numpy as _np
    dec0 = pd.read_csv(os.path.join(OUT, "fig4c_decay_vs_dose.csv"))
    dec1 = pd.read_csv(os.path.join(OUT, "fig4i_cycle_decay.csv"))
    con0 = pd.read_csv(os.path.join(OUT, "fig4d_fate_concordance.csv"))
    con1 = pd.read_csv(os.path.join(OUT, "fig4i_cycle_concordance.csv"))
    mel1 = pd.read_csv(os.path.join(OUT, "fig8f_melanoma_cycle.csv"))
    mel0 = pd.read_csv(os.path.join(OUT, "fig8a_melanoma_schedules.csv"))

    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.62, wspace=0.62)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(dec0.dose, dec0.decay_rate, marker="o", ms=3.5, lw=1.3, color=C[0], label="cycle-blind")
    ax.plot(dec1.dose, dec1.decay_rate, marker="s", ms=3.5, lw=1.3, color=C[1], label="cycle-gated")
    ax.set_xlabel("dose"); ax.set_ylabel("population decay rate")
    ax2 = ax.twinx(); ax2.grid(False)
    d0 = dec0[dec0.dose > 0]; d1 = dec1[dec1.dose > 0]
    ax2.plot(d0.dose, d0.death_time_mean, marker="o", ms=3, lw=1.0, ls="--", color=C[0], alpha=0.6)
    ax2.plot(d1.dose, d1.death_time_mean, marker="s", ms=3, lw=1.0, ls="--", color=C[1], alpha=0.6)
    ax2.set_ylabel("mean death time (dashed)", fontsize=6.5, labelpad=2); ax2.tick_params(labelsize=6.5)
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2)
    ax.set_title("a   decay and single-cell timing", loc="left", fontweight="semibold", color="#2a2a28")

    ax = fig.add_subplot(gs[0, 1])
    keys = [(m, r) for m in ("memory", "fast") for r in ("sisters", "cousins")]
    def excess(df, m, r):
        g = df[(df.model == m) & (df.relation == r)]
        return float(g.concordance.iloc[0] - g.expected_independent.iloc[0]) if len(g) else _np.nan
    x = _np.arange(len(keys)); w = 0.38
    ax.bar(x - w/2, [excess(con0, m, r) for m, r in keys], width=w, color=C[0], label="cycle-blind")
    ax.bar(x + w/2, [excess(con1, m, r) for m, r in keys], width=w, color=C[1], label="cycle-gated")
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"{m}\n{r}" for m, r in keys], fontsize=6)
    ax.set_ylabel("excess concordance")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2)
    ax.set_title("b   do related cells still share fates?", loc="left", fontweight="semibold", color="#2a2a28")

    ax = fig.add_subplot(gs[1, :])
    mechs = [m for m in dict.fromkeys(mel0.mechanism)]
    scheds = [s for s in dict.fromkeys(mel0.schedule)]
    x = _np.arange(len(mechs)); w = 0.8 / (2 * len(scheds))
    for j, s in enumerate(scheds):
        for k, (df, lab, alpha) in enumerate(((mel0, "cycle-blind", 1.0), (mel1, "cycle-gated", 0.45))):
            vals = [df[(df.mechanism == m) & (df.schedule == s)].ttp_baseline_weeks.mean() for m in mechs]
            sds = [df[(df.mechanism == m) & (df.schedule == s)].ttp_baseline_weeks.std() for m in mechs]
            off = (j * 2 + k - (2 * len(scheds) - 1) / 2) * w
            ax.bar(x + off, vals, yerr=sds, width=w, color=C[j], alpha=alpha, capsize=1.2,
                   label=(s if k == 0 else None))
    cens = mel0.ttp_baseline_weeks.max()
    ax.axhline(cens, color=INK2, ls=":", lw=0.8)
    ax.text(len(mechs) - 0.55, cens + 0.8, "end of follow-up", fontsize=5.5, color=INK2, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(["no fitness cost", "fitness cost", "partial protection"], fontsize=6.5)
    ax.set_ylabel("weeks to loss of control")
    h, l = ax.get_legend_handles_labels()
    h += [plt.Rectangle((0, 0), 1, 1, fc=INK2, alpha=a) for a in (1.0, 0.45)]; l += ["cycle-blind", "cycle-gated"]
    ax.legend(h, l, fontsize=6, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.16), columnspacing=1.0)
    ax.set_title("c   melanoma schedule ranking", loc="left", fontweight="semibold", color="#2a2a28")
    save(fig, "figS5_cycle_robustness")
except Exception as e:
    print("S5 skipped:", e)

# S6: the population and lineage layers against exact stationary laws
try:
    import numpy as _np
    val = pd.read_csv(os.path.join(OUT, "fig9_validation.csv"))
    pmf = pd.read_csv(os.path.join(OUT, "fig9_pmf.csv"))
    cnt = pd.read_csv(os.path.join(OUT, "fig9_counts.csv"))
    LABEL = {"constitutive": "constitutive", "bursty": "bursts of 4", "telegraph": "telegraph",
             "replication": "replication"}
    MODES = ("lineage", "population")
    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.42)

    # (a) simulated against exact distributions, both modes, for the constitutive model
    ax = fig.add_subplot(gs[0, 0])
    base = val[val.case == "constitutive"].dt.iloc[0]
    for i, mode in enumerate(MODES):
        g = cnt[(cnt.case == "constitutive") & (cnt["mode"] == mode) & (cnt.kernel == "DirectSSA") & (cnt.dt == base)]
        h = g.groupby("n").cells.sum()
        n = _np.arange(0, 60)
        freq = _np.array([h.get(k, 0) for k in n], dtype=float); freq /= freq.sum()
        ax.step(n, freq, where="mid", color=C[i], lw=1.1, label=f"{mode} (simulated)")
        e = pmf[pmf["mode"] == mode].set_index("n").exact
        ax.plot(n, [e.get(k, 0.0) for k in n], ls="--", lw=1.0, color=INK2,
                label="exact" if i == 0 else None)
    ax.set_xlabel("molecules per cell"); ax.set_ylabel("frequency")
    ax.legend(fontsize=6, loc="upper right")
    ax.set_title("a   one model, two modes", loc="left", fontweight="semibold", color="#2a2a28")

    # (b) distance to the exact law of the matching mode and of the other mode
    ax = fig.add_subplot(gs[0, 1])
    main = val[(val.case != "constitutive_dt") & (val.kernel == "DirectSSA")]
    keys = [(c, m) for c in ("constitutive", "bursty", "telegraph", "replication") for m in MODES
            if len(main[(main.case == c) & (main["mode"] == m)])]
    x = _np.arange(len(keys)); w = 0.38
    own = [main[(main.case == c) & (main["mode"] == m)].ks.max() for c, m in keys]
    oth = [main[(main.case == c) & (main["mode"] == m)].ks_other_mode.min() for c, m in keys]
    crit = [main[(main.case == c) & (main["mode"] == m)].ks_crit99.min() for c, m in keys]
    ax.bar(x - w / 2, own, width=w, color=C[0], label="to its own mode")
    ax.bar(x + w / 2, oth, width=w, color=C[1], label="to the other mode")
    for xi, ci in zip(x, crit):
        ax.plot([xi - 0.45, xi + 0.45], [ci, ci], ls=":", lw=0.9, color=INK2)
    ax.set_yscale("log"); ax.set_ylabel("Kolmogorov–Smirnov distance")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{LABEL[c]}\n{m[:4]}." for c, m in keys], fontsize=5.5)
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)
    ax.set_title("b   agreement, and its power to discriminate", loc="left", fontweight="semibold", color="#2a2a28")

    # (c) convergence in the update step
    ax = fig.add_subplot(gs[1, 0])
    dtv = val[val.case == "constitutive_dt"].sort_values("dt")
    ax.plot(dtv["dt"], dtv.ks_scheme_vs_continuum, marker="o", ms=3.5, lw=1.2, color=C[0],
            label="scheme against continuum (theory)")
    ax.plot(dtv["dt"], dtv.ks_continuum, marker="s", ms=3.5, lw=1.2, color=C[1],
            label="simulated against continuum")
    ax.plot(dtv["dt"], dtv.ks_crit99, ls=":", lw=0.9, color=INK2, label="99% critical value")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("update step / mean interdivision time"); ax.set_ylabel("Kolmogorov–Smirnov distance")
    ax.legend(fontsize=6)
    ax.set_title("c   the step-size bias is first order", loc="left", fontweight="semibold", color="#2a2a28")

    # (d) the emergent distribution of cell-cycle phase
    ax = fig.add_subplot(gs[1, 1])
    age = pd.read_csv(os.path.join(OUT, "fig9_agefit.csv"))
    for i, mode in enumerate(MODES):
        g = age[age["mode"] == mode].groupby("phase")[["observed", "expected"]].mean().reset_index()
        k = max(1, len(g) // 60)
        gg = g.iloc[::k]
        ax.plot(gg.phase, gg.observed / gg.observed.sum(), lw=1.1, color=C[i], label=f"{mode} (measured)")
        ax.plot(gg.phase, gg.expected / gg.expected.sum(), ls="--", lw=1.0, color=INK2,
                label="predicted" if i == 0 else None)
    ax.set_xlabel("cell-cycle phase"); ax.set_ylabel("fraction of cells")
    ax.legend(fontsize=6)
    ax.set_title("d   where the two modes differ", loc="left", fontweight="semibold", color="#2a2a28")
    save(fig, "figS6_exact")
except Exception as e:
    print("S6 skipped:", e)
