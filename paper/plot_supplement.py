"""Supplementary figures: S1 death/division-time distributions per dose, S2 ABC tolerance schedules, S3 runtime table."""
import os, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "output"); FIG = os.path.join(HERE, "figures")
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
plt.rcParams.update({"font.size": 7.5, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 300, "pdf.fonttype": 42})
def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight"); fig.savefig(os.path.join(FIG, name + ".png"), bbox_inches="tight"); plt.close(fig); print("saved", name)
# S1: single-cell timing distributions
try:
    d = pd.read_csv(os.path.join(OUT, "fig4c_times.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.4))
    for i, (dose, g) in enumerate(d[d.event == "division"].groupby("dose")):
        axes[0].hist(g.time, bins=np.linspace(10, 35, 40), histtype="step", color=C[i % 8], label=f"dose {dose:g}", density=True)
    axes[0].set_xlabel("division time (cells born under drug)"); axes[0].set_ylabel("density"); axes[0].legend(fontsize=6); axes[0].set_title("a  division times", loc="left")
    for i, (dose, g) in enumerate(d[d.event == "death"].groupby("dose")):
        axes[1].hist(g.time, bins=np.linspace(0, 60, 40), histtype="step", color=C[i % 8], label=f"dose {dose:g}", density=True)
    axes[1].set_xlabel("time from drug (or birth) to death"); axes[1].legend(fontsize=6); axes[1].set_title("b  death times", loc="left")
    save(fig, "figS1_timing")
except Exception as e: print("S1 skipped:", e)
# S2: ABC schedules
try:
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.4))
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
