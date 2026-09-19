"""Figure 1: schematic of the Biomodelling.jl 2.0 framework (drawn with matplotlib primitives)."""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Ellipse, Polygon
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, SURF = "#0b0b0b", "#52514e", "#f4f3ef"
plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42})
fig = plt.figure(figsize=(7.2, 5.9)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(4, 86); ax.axis("off")

def box(x, y, w, h, title, lines, color, title_size=8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=SURF, ec=color, lw=1.2))
    ax.text(x + 1.2, y + h - 1.4, title, fontsize=title_size, fontweight="bold", color=color, va="top")
    for i, l in enumerate(lines):
        ax.text(x + 1.2, y + h - 4.3 - 2.55 * i, l, fontsize=6.3, color=INK, va="top")
def arrow(x1, y1, x2, y2, color=INK2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9, lw=1.0, color=color))
def letter(x, y, s): ax.text(x, y, s, fontsize=10, fontweight="bold", va="top")

# (a) the cell
letter(1, 85.5, "a")
box(2, 58, 42, 26, "Cell: reactions in a growing, dividing volume", [], C[0])
# cell drawing: mother growing then dividing
ax.add_patch(Ellipse((10, 70), 9, 8, fc="white", ec=C[0], lw=1.2)); ax.text(10, 75.5, "V(t) = V₀ e^{λt}", ha="center", fontsize=6, color=INK2)
ax.text(10, 70.8, "G_off ⇄ G_on", ha="center", fontsize=5.6); ax.text(10, 68.6, "→ mRNA → protein", ha="center", fontsize=5.6)
ax.text(10, 66.6, "propensity ∝ V^(1−order)", ha="center", fontsize=5.2, color=INK2)
arrow(15, 70, 19.5, 70); ax.text(17.2, 71.2, "grow", ha="center", fontsize=5.5, color=INK2)
ax.add_patch(Ellipse((24, 70), 12, 9.5, fc="white", ec=C[0], lw=1.2)); ax.text(24, 70.6, "replication", ha="center", fontsize=5.6); ax.text(24, 68.4, "copies 1 → 2", ha="center", fontsize=5.6)
arrow(30.5, 70, 34.5, 72.5); arrow(30.5, 70, 34.5, 67.5); ax.text(32.3, 71.3, "divide", ha="center", fontsize=5.5, color=INK2, rotation=0)
ax.add_patch(Ellipse((38.5, 74), 8, 6.5, fc="white", ec=C[0], lw=1.2)); ax.add_patch(Ellipse((38.5, 66), 8, 6.5, fc="white", ec=C[0], lw=1.2))
ax.text(38.5, 74.4, "Binomial(n, f)", ha="center", fontsize=5.2); ax.text(38.5, 72.4, "promoters inherited", ha="center", fontsize=4.8, color=INK2)
ax.text(38.5, 66.4, "Binomial(n, 1−f)", ha="center", fontsize=5.2); ax.text(38.5, 64.4, "one copy each", ha="center", fontsize=4.8, color=INK2)
ax.text(3.2, 61.6, "kernels: DirectSSA · TauLeap · HybridSSATau · AdaptiveTauLeap        size control: Sizer · Adder · AgeTimer", fontsize=5.6, color=INK2)
ax.text(3.2, 59.4, "each cell: own RNG (deterministic, threaded) · parameters modulated per cell · volume held fixed within dt", fontsize=5.6, color=INK2)

# (b) population, drug, lineage
letter(47, 85.5, "b")
box(48, 58, 50, 26, "Population · drug · lineage", [
    "control: ConstantN (Moran replacement) · FreeGrowth · LogisticGrowth(K)",
    "dose d(t): constant · pulsed (holidays) · piecewise · bolus PK · function",
    "death hazard  h = h_max d^m/(EC50^m+d^m) · K^q/(K^q + c_P^q)",
    "growth inhibition · RateModulation(param, f(d)) · GenePerturbation",
    "LineageTable: parent · clone · generation · birth/end time ·",
    "    fate · state at birth and at division → Newick",
], C[1])
# tiny tree
tx, ty = 90, 60.0
for (x1, y1, x2, y2) in ((0, 0, -4, 3), (0, 0, 4, 3), (-4, 3, -6, 6), (-4, 3, -2, 6), (4, 3, 2, 6), (4, 3, 6, 6)):
    ax.plot([tx + x1, tx + x2], [ty + y1, ty + y2], color=C[1], lw=1.0)
for (x, y, c) in ((-6, 6, C[7]), (-2, 6, C[1]), (2, 6, C[1]), (6, 6, C[7])): ax.add_patch(Circle((tx + x, ty + y), 0.9, fc=c, ec="none"))
ax.text(tx, ty - 2.2, "fates: ● divided  ● died", ha="center", fontsize=5, color=INK2)

# (c) observation
letter(1, 55, "c")
box(2, 33, 42, 20, "Observation and output", [
    "scRNA-seq: capture ~ Beta(mean, cv) · Poisson depth ·",
    "   dropout · batches; library size ∝ cell volume",
    "smFISH: fixed detection efficiency",
    "time-lapse: reporter sampled along lineages",
    "outputs: CSV · Newick · AnnData (.h5ad) with true",
    "   volume, age, clone and copy number per cell",
], C[2])

# (d) downstream
letter(47, 55, "d")
box(48, 33, 50, 20, "Downstream", [
    "lineage statistics: mother–daughter, sister and cousin correlations,",
    "   lineage autocorrelation, memory timescale, population vs lineage noise",
    "fluctuation tests (Luria–Delbrück) · clonal (MemorySeq-style) scores",
    "inference: ABC-SMC over any simulation · exact telegraph likelihood",
    "generators: random GRNs (ER / scale-free) with known ground truth",
], C[6])
arrow(23, 57.5, 23, 53.5); arrow(73, 57.5, 73, 53.5); arrow(44.5, 43, 47.5, 43)

# (e) code
letter(1, 30, "e")
code = """using Biomodelling, Random
model = telegraph_model(k_on = 0.005, k_off = 0.005, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)   # resistance gene R → P
x0    = initial_state(model; G_off = 1)
pert  = Perturbation(PiecewiseDose([0.0, 100.0], [0.0, 1.0]);                                          # drug from t = 100
                     effects = [DeathHazard(h_max = 0.5, EC50 = 0.3, protect = :protein, K = 60.0, q = 4.0)])
st    = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2)/20), size_control = Sizer(2.0; cv = 0.05),
                           replication = Replication(0.5), control = FreeGrowth(max_cells = 20_000))
res   = simulate_population(model, x0, 1000, (0.0, 200.0); settings = st, perturbation = pert, rng = Xoshiro(1))
heritability(res, :protein; relation = :mother_daughter)              # memory across divisions
Y     = sequence(sample_cells(res; n = 500).counts, SeqProtocol(capture = 0.15)).Y   # synthetic scRNA-seq"""
ax.add_patch(FancyBboxPatch((2, 6), 96, 22, boxstyle="round,pad=0.4,rounding_size=1.2", fc="#fbfbfa", ec=INK2, lw=0.8))
ax.text(3.2, 26.8, code, fontsize=5.6, family="monospace", va="top", color=INK, linespacing=1.35)
fig.savefig(os.path.join(FIG, "fig1_framework.pdf"), bbox_inches="tight"); fig.savefig(os.path.join(FIG, "fig1_framework.png"), bbox_inches="tight", dpi=300)
print("saved fig1_framework")
