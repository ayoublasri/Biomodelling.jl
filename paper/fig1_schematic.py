"""Figure 1: schematic of the Biomodelling.jl 2.0 framework (vector drawing with matplotlib primitives)."""
import os, textwrap, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Ellipse, Rectangle
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, SURF = "#0b0b0b", "#52514e", "#f6f5f2"
plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "mathtext.fontset": "dejavusans"})
W, H = 100.0, 90.0
fig = plt.figure(figsize=(7.2, 5.9)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(10, H); ax.axis("off")
def fx(x): return x / W
def fy(y): return (y - 10) / (H - 10)

def box(x, y, w, h, title, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=SURF, ec=color, lw=1.2))
    ax.text(x + 1.3, y + h - 1.3, title, fontsize=7.6, fontweight="bold", color=color, va="top")
def lines(x, y, items, dy=2.5, size=6.2, color=INK):
    for i, l in enumerate(items): ax.text(x, y - dy * i, l, fontsize=size, color=color, va="top")
def arrow(x1, y1, x2, y2, color=INK2, lw=1.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9, lw=lw, color=color))
def letter(x, y, s): ax.text(x, y, s, fontsize=10, fontweight="bold", va="top")
def bg(txt): return dict(boxstyle="round,pad=0.15", fc=SURF, ec="none")

# ---------------------------------------------------------------- (a) the cell
letter(1, 89.3, "a")
box(2, 61, 43, 27, "Cell: reactions in a growing, dividing volume", C[0])
ax.add_patch(Ellipse((10.5, 75), 10, 8.6, fc="white", ec=C[0], lw=1.2))
ax.text(10.5, 76.4, r"$G_{\mathrm{off}} \rightleftharpoons G_{\mathrm{on}}$", ha="center", fontsize=6.4)
ax.text(10.5, 73.6, "mRNA → protein", ha="center", fontsize=5.8)
ax.text(10.5, 81.2, r"$V(t) = V_0\,e^{\lambda t}$", ha="center", fontsize=6.2, color=INK2)
ax.text(10.5, 69.2, r"propensity $\propto V^{\,1-\mathrm{order}}$", ha="center", fontsize=5.4, color=INK2)
arrow(16.0, 75, 19.6, 75); ax.text(17.8, 76.1, "grow", ha="center", fontsize=5.4, color=INK2)
ax.add_patch(Ellipse((26, 75), 12, 10.2, fc="white", ec=C[0], lw=1.2))
ax.text(26, 76.3, "gene replication", ha="center", fontsize=5.8); ax.text(26, 73.8, "copies 1 → 2", ha="center", fontsize=5.8)
arrow(32.4, 75.6, 35.4, 78.6); arrow(32.4, 74.4, 35.4, 71.4); ax.text(34.6, 75.0, "divide", ha="center", va="center", fontsize=5.2, color=INK2)
ax.add_patch(Ellipse((40, 80.2), 8.4, 6.2, fc="white", ec=C[0], lw=1.2)); ax.add_patch(Ellipse((40, 69.8), 8.4, 6.2, fc="white", ec=C[0], lw=1.2))
ax.text(40, 80.9, "Binomial(n, f)", ha="center", fontsize=5.3); ax.text(40, 78.8, "promoters inherited", ha="center", fontsize=4.6, color=INK2)
ax.text(40, 70.5, "Binomial(n, 1−f)", ha="center", fontsize=5.3); ax.text(40, 68.4, "one gene copy each", ha="center", fontsize=4.6, color=INK2)
lines(3.3, 65.6, ["kernels: DirectSSA · TauLeap · HybridSSATau · AdaptiveTauLeap (volume fixed within dt)",
                  "size control: Sizer · Adder · AgeTimer   ·   per-cell RNG and parameters"], dy=2.3, size=5.5, color=INK2)

# ---------------------------------------------------------------- (b) population, drug, lineage
letter(47, 89.3, "b")
box(48, 61, 50, 27, "Population, drug and lineage", C[1])
lines(49.3, 84.0, ["control: constant N (Moran replacement) · free growth · logistic",
                   "dose d(t): constant · pulsed · piecewise · bolus PK · function",
                   "effects: death hazard · growth inhibition · rate modulation · fitness cost · gene perturbation"], dy=2.3, size=5.5)
# dose schedule icon
axd = fig.add_axes([fx(50.5), fy(64.0), 0.13, 0.10]); t = np.linspace(0, 100, 400); d = ((t > 30) & (((t - 30) % 30) < 20)).astype(float)
axd.plot(t, d, color=C[1], lw=1.2); axd.set_ylim(-0.1, 1.5); axd.set_xticks([]); axd.set_yticks([]); axd.set_facecolor("white")
for s in axd.spines.values(): s.set_color(INK2); s.set_linewidth(0.6)
axd.set_title("dose schedule d(t)", fontsize=5.2, color=INK2, pad=2)
axd.set_xlabel("time", fontsize=5, labelpad=1)
# hazard icon
axh = fig.add_axes([fx(67.5), fy(64.0), 0.11, 0.10]); dd = np.linspace(0, 2, 200); h = dd**2 / (0.25 + dd**2)
axh.plot(dd, h, color=C[7], lw=1.2, label="unprotected"); axh.plot(dd, h / 20, color=C[0], lw=1.2, label="protected")
axh.set_xticks([]); axh.set_yticks([]); axh.set_facecolor("white"); axh.legend(fontsize=4.4, frameon=False, loc="upper left", handlelength=1.2)
for s in axh.spines.values(): s.set_color(INK2); s.set_linewidth(0.6)
axh.set_title(r"death hazard $h(d, c_P)$", fontsize=5.2, color=INK2, pad=2); axh.set_xlabel("dose", fontsize=5, labelpad=1)
# lineage tree icon
tx, ty = 89.5, 66.5
for (x1, y1, x2, y2) in ((0, 0, -4, 3), (0, 0, 4, 3), (-4, 3, -6.2, 6.5), (-4, 3, -1.8, 6.5), (4, 3, 1.8, 6.5), (4, 3, 6.2, 6.5)):
    ax.plot([tx + x1, tx + x2], [ty + y1, ty + y2], color=C[1], lw=1.0, solid_capstyle="round")
for (x, y, c) in ((-6.2, 6.5, C[7]), (-1.8, 6.5, C[1]), (1.8, 6.5, C[1]), (6.2, 6.5, C[7])): ax.add_patch(Circle((tx + x, ty + y), 0.85, fc=c, ec="none"))
ax.text(tx, ty + 8.4, "lineage table", ha="center", fontsize=5.2, color=INK2)
ax.text(tx, ty - 1.6, "parent · clone · birth/end time", ha="center", fontsize=4.6, color=INK2)
ax.text(tx, ty - 3.4, "fate (● divided ● died) · state", ha="center", fontsize=4.6, color=INK2)

# ---------------------------------------------------------------- (c) observation
letter(1, 58.3, "c")
box(2, 35, 43, 22, "Observation and output", C[2])
# count matrix icon: true counts -> sequenced counts
rng = np.random.default_rng(0)
def matrix(x0, y0, w, h, vals, cmap):
    n, m = vals.shape
    for i in range(n):
        for j in range(m):
            ax.add_patch(Rectangle((x0 + j * w / m, y0 + i * h / n), w / m, h / n, fc=cmap(vals[i, j]), ec="white", lw=0.3))
true = rng.gamma(2.0, 1.0, (6, 8)); true /= true.max()
seq = np.where(rng.random((6, 8)) < 0.35, 0, true * rng.beta(2, 8, (6, 8)) * 4); seq = np.clip(seq / max(seq.max(), 1e-9), 0, 1)
matrix(4, 41.5, 10, 8, true, plt.get_cmap("Blues")); ax.text(9, 50.3, "true counts", ha="center", fontsize=5.2, color=INK2)
ax.text(9, 40.5, "cells × genes", ha="center", fontsize=4.6, color=INK2)
arrow(15, 45.5, 19.5, 45.5); ax.text(17.2, 46.6, "observe", ha="center", fontsize=5, color=INK2)
matrix(20.5, 41.5, 10, 8, seq, plt.get_cmap("Greens")); ax.text(25.5, 50.3, "scRNA-seq counts", ha="center", fontsize=5.2, color=INK2)
ax.text(25.5, 40.5, "capture ~ Beta · depth · dropout · batch", ha="center", fontsize=4.4, color=INK2)
lines(32.5, 49.6, ["smFISH: detection", "   efficiency", "time-lapse: reporter", "   along lineages"], dy=2.1, size=5.4)
lines(3.3, 38.6, ["outputs: CSV · Newick · AnnData (.h5ad) with volume, age, clone and copy number"], size=5.5, color=INK2)

# ---------------------------------------------------------------- (d) downstream
letter(47, 58.3, "d")
box(48, 35, 50, 22, "Downstream analyses", C[6])
items = [("lineage statistics", "mother–daughter, sister and cousin correlations; lineage autocorrelation and memory timescale; population vs lineage noise"),
         ("clonal tests", "Luria–Delbrück fluctuation test; MemorySeq-style clonal variance scores"),
         ("inference", "ABC-SMC with any simulation as forward model; exact telegraph likelihood"),
         ("benchmarks", "random GRNs (Erdős–Rényi, scale-free) with known ground truth; knockdown / overexpression datasets")]
y = 53.4
for i, (k, v) in enumerate(items):
    ax.add_patch(Circle((50.2, y - 0.6), 0.55, fc=C[6], ec="none"))
    ax.text(51.4, y, k, fontsize=5.9, fontweight="bold", color=C[6], va="top")
    ax.text(51.4, y - 2.2, "\n".join(textwrap.wrap(v, 86)), fontsize=5.2, color=INK, va="top", linespacing=1.15)
    y -= 4.6
arrow(23.5, 60.5, 23.5, 57.4); arrow(73, 60.5, 73, 57.4); arrow(45.5, 46, 47.5, 46)

# ---------------------------------------------------------------- (e) code
letter(1, 32.3, "e")
code = """using Biomodelling, Random
model = telegraph_model(k_on = 0.002, k_off = 0.02, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)   # resistance gene R → mRNA → P
x0    = initial_state(model; G_off = 1)
pert  = Perturbation(PiecewiseDose([0.0, 40.0], [0.0, 1.0]);                                          # drug from t = 40
                     effects = [DeathHazard(h_max = 0.15, EC50 = 0.5, m = 2.0, protect = :protein, K = 150.0, q = 4.0)])
st    = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2)/20), size_control = Sizer(2.0; cv = 0.05),
                           replication = Replication(0.5), control = FreeGrowth(max_cells = 20_000))
res   = simulate_population(model, x0, 500, (0.0, 200.0); settings = st, perturbation = pert, rng = Xoshiro(1))
heritability(res, :protein; relation = :sisters)                       # expression memory across divisions
Y     = sequence(sample_cells(res; n = 500).counts, SeqProtocol(capture = 0.15)).Y    # synthetic scRNA-seq"""
ax.add_patch(FancyBboxPatch((2, 12), 96, 17.5, boxstyle="round,pad=0.4,rounding_size=1.2", fc="#fbfbfa", ec=INK2, lw=0.8))
ax.text(4.6, 30.4, "A complete example: a resistance gene in a growing, dividing, drug-treated population", fontsize=6.2, fontweight="bold", color=INK2, va="bottom")
ax.text(3.2, 28.6, code, fontsize=5.6, family="monospace", va="top", color=INK, linespacing=1.4)
fig.savefig(os.path.join(FIG, "fig1_framework.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(FIG, "fig1_framework.png"), bbox_inches="tight", dpi=300)
print("saved fig1_framework")
