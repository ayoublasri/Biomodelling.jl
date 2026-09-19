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
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
    for i, (dose, g) in enumerate(d[d.event == "division"].groupby("dose")):
        axes[0].hist(g.time, bins=np.linspace(10, 35, 40), histtype="step", color=C[i % 8], label=f"dose {dose:g}", density=True)
    axes[0].set_xlabel("division time (cells born under drug)"); axes[0].set_ylabel("density"); axes[0].legend(fontsize=6); axes[0].set_title("a   division times", loc="left", fontweight="semibold", color="#2a2a28")
    for i, (dose, g) in enumerate(d[d.event == "death"].groupby("dose")):
        axes[1].hist(g.time, bins=np.linspace(0, 60, 40), histtype="step", color=C[i % 8], label=f"dose {dose:g}", density=True)
    axes[1].set_xlabel("time from drug (or birth) to death"); axes[1].legend(fontsize=6); axes[1].set_title("b   death times", loc="left", fontweight="semibold", color="#2a2a28")
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
