"""Figure 1: the Biomodelling.jl 2.0 framework.

Drawn as an SVG document and rendered with headless Chromium, which gives control
over typography and spacing that the matplotlib primitives did not:

    python fig1_schematic.py                       # writes figures/fig1_framework.html
    chromium --headless --disable-gpu --hide-scrollbars --force-device-scale-factor=3 \
        --window-size=1440,1060 --virtual-time-budget=4000 \
        --screenshot=figures/fig1_raw.png figures/fig1_framework.html
    # crop the top 4320x2940 pixels to figures/fig1_framework.png (600 dpi at 7.2 in)

`make_figure1()` returns the HTML; running the module writes it next to the figures.
"""
import base64, os, re, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
FONTS = {  # Inter, served by Google Fonts; cached in figures/fonts/ on first use
    "400": "https://fonts.gstatic.com/s/inter/v20/UcCO3FwrK3iLTeHuS_nVMrMxCp50SjIw2boKoduKmMEVuLyfMZg.ttf",
    "600": "https://fonts.gstatic.com/s/inter/v20/UcCO3FwrK3iLTeHuS_nVMrMxCp50SjIw2boKoduKmMEVuGKYMZg.ttf",
    "800": "https://fonts.gstatic.com/s/inter/v20/UcCO3FwrK3iLTeHuS_nVMrMxCp50SjIw2boKoduKmMEVuDyYMZg.ttf",
}

# palette shared with the data figures
BLUE, ORANGE, GREEN, AMBER, PINK, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#e34948"
INK, INK2, MUTED, LINE, PANEL, PAPER = "#0b0b0b", "#52514e", "#8a8985", "#e6e5e1", "#f7f6f3", "#ffffff"
TYPE = 1.52   # every font size is multiplied by this so nothing prints below ~5 pt


def font(weight):
    path = os.path.join(FIG, "fonts", f"Inter-{weight}.ttf")
    if not os.path.exists(path):
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            urllib.request.urlretrieve(FONTS[str(weight)], path)
        except Exception as exc:                                  # offline: fall back to a system sans
            print(f"  [Inter {weight} unavailable ({exc}); using the system sans-serif]")
            return ""
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


# ---------------------------------------------------------------- drawing helpers
def cell(cx, cy, r, state, *, mrna=0, protein=0, copies=1, dead=False, seed=0):
    """A cell: membrane, nucleus with 1 or 2 gene copies, mRNA and protein dots."""
    import random
    rng = random.Random(seed)
    edge = MUTED if dead else (ORANGE if state == "on" else BLUE)
    body = "#fdfdfc" if not dead else "#f0efec"
    dash = ' stroke-dasharray="3 3"' if dead else ''
    out = [f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{body}" stroke="{edge}" stroke-width="1.6"{dash}/>']
    for i in range(mrna):                                          # mRNA: small open squares
        a, d = rng.uniform(0, 6.283), rng.uniform(0.42, 0.86) * r
        x, y = cx + d * __import__("math").cos(a), cy + d * __import__("math").sin(a)
        out.append(f'<rect x="{x-2:.1f}" y="{y-2:.1f}" width="4" height="4" rx="1" fill="{GREEN}"/>')
    for i in range(protein):                                       # protein: filled dots
        a, d = rng.uniform(0, 6.283), rng.uniform(0.45, 0.9) * r
        x, y = cx + d * __import__("math").cos(a), cy + d * __import__("math").sin(a)
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.6" fill="{VIOLET}"/>')
    out.append(f'<circle cx="{cx}" cy="{cy}" r="{r*0.42:.1f}" fill="#ffffff" stroke="{LINE}" stroke-width="1"/>')
    gx = [cx] if copies == 1 else [cx - r * 0.16, cx + r * 0.16]    # gene copies in the nucleus
    for x in gx:
        out.append(f'<rect x="{x-3.4:.1f}" y="{cy-3.4:.1f}" width="6.8" height="6.8" rx="1.6" '
                   f'fill="{ORANGE if state == "on" else "#ffffff"}" stroke="{ORANGE if state == "on" else MUTED}" stroke-width="1.4"/>')
    if dead:
        out.append(f'<path d="M{cx-r*0.5:.1f} {cy-r*0.5:.1f} L{cx+r*0.5:.1f} {cy+r*0.5:.1f} '
                   f'M{cx+r*0.5:.1f} {cy-r*0.5:.1f} L{cx-r*0.5:.1f} {cy+r*0.5:.1f}" stroke="{RED}" stroke-width="2" stroke-linecap="round"/>')
    return "".join(out)


def arrow(x1, y1, x2, y2, colour=INK2, width=1.4, dash=None, head=True):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    marker = ' marker-end="url(#arrowhead)"' if head else ""
    return f'<path d="M{x1} {y1} L{x2} {y2}" stroke="{colour}" stroke-width="{width}" fill="none"{d}{marker}/>'


def curve(x1, y1, x2, y2, bend=30, colour=INK2, width=1.4, head=True):
    marker = ' marker-end="url(#arrowhead)"' if head else ""
    return (f'<path d="M{x1} {y1} Q{(x1+x2)/2} {(y1+y2)/2 - bend} {x2} {y2}" '
            f'stroke="{colour}" stroke-width="{width}" fill="none"{marker}/>')


def label(x, y, text, size=12.5, weight=400, colour=INK, anchor="start", cls=""):
    return (f'<text x="{x}" y="{y}" font-size="{size*TYPE:.1f}" font-weight="{weight}" fill="{colour}" '
            f'text-anchor="{anchor}" class="{cls}">{text}</text>')


def panel(x, y, w, h, letter, title, accent):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{PANEL}" stroke="{LINE}" stroke-width="1"/>'
            f'<rect x="{x}" y="{y}" width="4" height="{h}" rx="2" fill="{accent}"/>'
            + label(x + 18, y + 30, letter, 15.5, 800, INK)
            + label(x + 44, y + 30, title, 11.6, 600, INK))


def make_figure1():
    W, H = 1440, 1060
    s = []

    # ---------------------------------------------------------------- a: one cell
    s.append(panel(24, 20, 806, 300, "a", "One cell: stochastic kinetics in a growing, dividing volume", BLUE))
    s.append(label(60, 68, "Propensities scale with volume, so concentrations are invariant under growth.", 11.5, 400, INK2))

    s.append(cell(104, 172, 31, "off", mrna=2, protein=3, seed=1))
    s.append(label(104, 228, "birth", 11, 600, INK, "middle"))
    s.append(label(104, 248, "one gene copy", 10, 400, MUTED, "middle"))
    s.append(arrow(142, 172, 186, 172))
    s.append(label(164, 158, "grow", 10.5, 600, INK2, "middle"))

    s.append(cell(244, 172, 42, "on", mrna=5, protein=7, seed=2))
    s.append(label(244, 228, "gene replication", 11, 600, INK, "middle"))
    s.append(label(244, 248, "copies 1 → 2", 10, 400, MUTED, "middle"))
    s.append(arrow(292, 172, 338, 172))
    s.append(label(315, 158, "divide", 10.5, 600, INK2, "middle"))

    s.append(cell(392, 132, 28, "on", mrna=3, protein=4, seed=3))
    s.append(cell(392, 212, 28, "off", mrna=2, protein=3, seed=4))
    s.append(label(432, 126, "molecules split Binomial(n, f)", 10, 400, INK))
    s.append(label(432, 145, "promoter state inherited", 10, 400, INK))
    s.append(label(432, 206, "volume split f ~ 𝒩(0.5, σ²)", 10, 400, INK))
    s.append(label(432, 225, "sizer · adder · timer rules", 10, 400, INK))

    # growth curve inset
    s.append(f'<rect x="654" y="102" width="154" height="120" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    s.append(label(666, 121, "volume", 10.5, 600, INK2))
    pts, x0, y0 = [], 668, 200
    for i in range(61):
        t = i / 60
        seg = t * 3 % 1
        pts.append(f"{x0 + t*126:.1f},{y0 - 32*(2**seg - 1) - 6:.1f}")
    s.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{BLUE}" stroke-width="1.8"/>')
    for k in (1, 2):
        xx = x0 + 126 * k / 3
        s.append(f'<path d="M{xx:.1f} 164 L{xx:.1f} 202" stroke="{MUTED}" stroke-width="1" stroke-dasharray="3 3"/>')
    s.append(label(668, 218, "halved at division", 10, 400, INK2))

    # legend of species
    s.append(f'<rect x="60" y="272" width="748" height="32" rx="7" fill="{PAPER}" stroke="{LINE}"/>')
    items = [("rect", ORANGE, "promoter (on / off)"), ("sq", GREEN, "mRNA"),
             ("dot", VIOLET, "protein"), ("x", RED, "killed by the drug")]
    lx = 76
    for kind, colour, text in items:
        if kind == "rect":
            s.append(f'<rect x="{lx}" y="284" width="9" height="9" rx="2" fill="{colour}"/>')
        elif kind == "sq":
            s.append(f'<rect x="{lx}" y="284" width="8" height="8" rx="1.5" fill="{colour}"/>')
        elif kind == "dot":
            s.append(f'<circle cx="{lx+4}" cy="288.5" r="4.4" fill="{colour}"/>')
        else:
            s.append(f'<path d="M{lx} 284 L{lx+9} 293 M{lx+9} 284 L{lx} 293" stroke="{colour}" stroke-width="1.8" stroke-linecap="round"/>')
        s.append(label(lx + 17, 293, text, 10, 400, INK))
        lx += 42 + int(len(text) * 8.3)
    return W, H, s


def add_panel_b(s):
    """b: population, drug schedule and the lineage record."""
    import math, random
    s.append(panel(854, 20, 562, 300, "b", "Population under drug, with its lineage", ORANGE))
    s.append(label(890, 68, "Dose schedules and state-dependent killing, applied to every", 11.5, 400, INK2))
    s.append(label(890, 85, "cell at every step, on a free, fixed or logistic population.", 11.5, 400, INK2))

    # dose schedule
    x0, y0, w = 890, 142, 214
    s.append(label(x0, 112, "dose d(t)", 11, 600, INK))
    steps = [(0, 0), (0.18, 0), (0.18, 1), (0.42, 1), (0.42, 0), (0.60, 0), (0.60, 1), (0.84, 1), (0.84, 0), (1, 0)]
    pts = " ".join(f"{x0 + p*w:.1f},{y0 - 26*v:.1f}" for p, v in steps)
    s.append(f'<polyline points="{pts}" fill="none" stroke="{ORANGE}" stroke-width="2"/>')
    for a, b in ((0.18, 0.42), (0.60, 0.84)):
        s.append(f'<rect x="{x0+a*w:.1f}" y="{y0-26:.1f}" width="{(b-a)*w:.1f}" height="26" fill="{ORANGE}" opacity="0.13"/>')
    s.append(f'<path d="M{x0} {y0} L{x0+w} {y0}" stroke="{LINE}" stroke-width="1"/>')
    s.append(label(x0, y0 + 20, "pulsed · piecewise · bolus PK", 10, 400, MUTED))

    # hazard sketch
    hx = 1148
    s.append(label(hx, 112, "death hazard", 11, 600, INK))
    s.append(f'<rect x="{hx}" y="{122}" width="240" height="46" rx="7" fill="{PAPER}" stroke="{LINE}"/>')
    s.append(f'<text x="{hx+12}" y="144" font-size="14" fill="{INK}">h = h<tspan font-size="10" dy="2">max</tspan>'
             f'<tspan font-size="14" dy="-2"> · f(dose) · g(protein)</tspan></text>')
    s.append(label(hx + 12, 162, "growth inhibition · fitness cost", 10, 400, MUTED))

    # lineage tree with fates
    tx, ty, tw, th = 890, 208, 420, 80
    s.append(label(tx, 198, "lineage record: ancestry, state at birth, fate", 11, 600, INK))
    rng = random.Random(7)
    gens = 4
    nodes = {0: [(tx + 10, ty + th / 2, "off")]}
    edges = []
    for g in range(1, gens):
        parents, children = nodes[g - 1], []
        span = th / (2 ** g)
        for (px, py, st) in parents:
            for k in (-1, 1):
                cy = py + k * span / 2
                cx = tx + 10 + g * (tw - 40) / (gens - 1)
                st2 = st if rng.random() > 0.22 else ("on" if st == "off" else "off")
                children.append((cx, cy, st2))
                edges.append((px, py, cx, cy))
        nodes[g] = children
    for (px, py, cx, cy) in edges:
        s.append(f'<path d="M{px:.1f} {py:.1f} L{px+14:.1f} {py:.1f} L{px+14:.1f} {cy:.1f} L{cx:.1f} {cy:.1f}" '
                 f'fill="none" stroke="{LINE if False else "#cfcdc7"}" stroke-width="1.3"/>')
    for g, pts in nodes.items():
        for i, (x, y, st) in enumerate(pts):
            killed = g == gens - 1 and st == "off" and rng.random() < 0.72
            colour = ORANGE if st == "on" else BLUE
            if killed:
                s.append(f'<path d="M{x-4:.1f} {y-4:.1f} L{x+4:.1f} {y+4:.1f} M{x+4:.1f} {y-4:.1f} L{x-4:.1f} {y+4:.1f}" '
                         f'stroke="{RED}" stroke-width="1.8" stroke-linecap="round"/>')
            else:
                s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.6" fill="{colour}"/>')
    s.append(f'<rect x="{tx+tw+8}" y="{ty-6}" width="1" height="{th+12}" fill="{LINE}"/>')
    s.append(label(tx + tw + 18, ty + 36, "surviving", 10, 400, INK2))
    s.append(label(tx + tw + 18, ty + 50, "clones", 10, 400, INK2))
    s.append(label(tx, ty + th + 24, "heritability · memory time · clone diversity · fluctuation tests", 10, 400, MUTED))


def add_panel_c(s):
    """c: observation models."""
    s.append(panel(24, 336, 620, 244, "c", "From molecules to measurements", GREEN))
    s.append(label(60, 386, "One true state, three kinds of measurement.", 11.5, 400, INK2))

    import random
    rng = random.Random(3)
    # true counts matrix -> sequenced matrix
    def matrix(x, y, cols, rows, cell_w, dropout):
        out = []
        for i in range(rows):
            for j in range(cols):
                v = rng.random()
                if dropout and rng.random() < 0.45:
                    fill, op = "#ffffff", 1.0
                else:
                    fill, op = GREEN, 0.18 + 0.8 * v
                out.append(f'<rect x="{x + j*cell_w:.1f}" y="{y + i*cell_w:.1f}" width="{cell_w-1.4:.1f}" '
                           f'height="{cell_w-1.4:.1f}" rx="1.5" fill="{fill}" opacity="{op:.2f}" stroke="{LINE}" stroke-width="0.5"/>')
        return "".join(out)

    s.append(label(60, 416, "true counts", 10.5, 600, INK))
    s.append(matrix(60, 426, 8, 7, 12, False))
    s.append(arrow(160, 468, 196, 468))
    s.append(label(178, 456, "capture", 10, 600, INK2, "middle"))

    s.append(label(208, 416, "scRNA-seq", 10.5, 600, INK))
    s.append(matrix(208, 426, 8, 7, 12, True))

    # smFISH
    s.append(label(340, 416, "smFISH", 10.5, 600, INK))
    s.append(f'<rect x="340" y="426" width="96" height="92" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    s.append(f'<circle cx="388" cy="472" r="33" fill="#fdfdfc" stroke="{BLUE}" stroke-width="1.4"/>')
    for _ in range(14):
        import math
        a, d = rng.uniform(0, 6.283), rng.uniform(0.15, 0.85) * 30
        s.append(f'<circle cx="{388 + d*math.cos(a):.1f}" cy="{472 + d*math.sin(a):.1f}" r="2.4" fill="{GREEN}"/>')

    # time-lapse reporter
    s.append(label(452, 416, "time-lapse", 10.5, 600, INK))
    s.append(f'<rect x="452" y="426" width="148" height="92" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    px, py = 462, 506
    pts = []
    val = 18
    for i in range(35):
        val += rng.uniform(-4, 4.4)
        val = max(4, min(62, val))
        pts.append(f"{px + i*3.75:.1f},{py - val*0.94:.1f}")
    s.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{VIOLET}" stroke-width="1.5"/>')
    for k in (1, 2):
        xx = px + 35 * 3.75 * k / 3
        s.append(f'<path d="M{xx:.1f} 436 L{xx:.1f} 508" stroke="{MUTED}" stroke-width="0.9" stroke-dasharray="3 3"/>')
    s.append(label(526, 534, "divisions", 10, 400, MUTED, "middle"))

    s.append(f'<rect x="60" y="542" width="560" height="30" rx="7" fill="{PAPER}" stroke="{LINE}"/>')
    s.append(label(74, 562, "output: CSV · AnnData · Newick · lineage tables", 10, 400, INK))


def add_panel_d(s):
    """d: what the framework is used for."""
    import random, math
    rng = random.Random(11)
    s.append(panel(668, 336, 748, 244, "d", "What the forward model is used for", VIOLET))

    # d1 memory genes
    x, y = 700, 394
    s.append(label(x, y, "memory genes", 10.5, 600, INK))
    s.append(f'<rect x="{x}" y="{y+10}" width="140" height="92" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    for i in range(26):
        t = i / 25
        sx = x + 13 + t * 114
        sy = y + 94 - (6 + 72 / (1 + math.exp(-(t - 0.55) * 11))) * 0.82
        s.append(f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="2.4" fill="{GREEN}" opacity="0.85"/>')
    s.append(f'<path d="M{x+75:.1f} {y+16} L{x+75:.1f} {y+96}" stroke="{MUTED}" stroke-width="1" stroke-dasharray="3 3"/>')
    s.append(label(x + 70, y + 120, "memory vs cycle", 10, 400, MUTED, "middle"))

    # d2 network inference
    x = 858
    s.append(label(x, y, "network inference", 10.5, 600, INK))
    s.append(f'<rect x="{x}" y="{y+10}" width="140" height="92" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    nodes = [(x + 38, y + 32), (x + 92, y + 28), (x + 112, y + 62), (x + 70, y + 82), (x + 28, y + 66)]
    for i, (nx, ny) in enumerate(nodes):
        for j in (1, 3):
            tx2, ty2 = nodes[(i + j) % len(nodes)]
            s.append(f'<path d="M{nx} {ny} L{tx2} {ty2}" stroke="#cfcdc7" stroke-width="1.1"/>')
    for i, (nx, ny) in enumerate(nodes):
        s.append(f'<circle cx="{nx}" cy="{ny}" r="6.4" fill="{[BLUE, ORANGE, GREEN, AMBER, VIOLET][i]}"/>')
    s.append(label(x + 70, y + 120, "known edges", 10, 400, MUTED, "middle"))

    # d3 ABC posterior
    x = 1016
    s.append(label(x, y, "parameters", 10.5, 600, INK))
    s.append(f'<rect x="{x}" y="{y+10}" width="140" height="92" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    for _ in range(110):
        u, v = rng.gauss(0, 1), rng.gauss(0, 1)
        px = x + 70 + 24 * u + 12 * v
        py = y + 58 + 16 * v
        if x + 8 < px < x + 132 and y + 16 < py < y + 96:
            s.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="2.1" fill="{BLUE}" opacity="0.5"/>')
    s.append(f'<path d="M{x+63} {y+50} L{x+77} {y+64} M{x+77} {y+50} L{x+63} {y+64}" stroke="{INK}" stroke-width="2" stroke-linecap="round"/>')
    s.append(label(x + 70, y + 120, "posterior, truth", 10, 400, MUTED, "middle"))

    # d4 schedule optimisation
    x = 1174
    s.append(label(x, y, "schedules", 10.5, 600, INK))
    s.append(f'<rect x="{x}" y="{y+10}" width="208" height="92" rx="8" fill="{PAPER}" stroke="{LINE}"/>')
    cw, ch = 25, 19
    best = (3, 1)
    for i in range(4):
        for j in range(7):
            d2 = ((j - best[0]) / 3.2) ** 2 + ((i - best[1]) / 1.9) ** 2
            v = max(0.06, 1 - d2 * 0.55)
            s.append(f'<rect x="{x + 12 + j*cw:.1f}" y="{y + 20 + i*ch:.1f}" width="{cw-2:.1f}" height="{ch-2:.1f}" '
                     f'rx="2" fill="{VIOLET}" opacity="{0.10 + 0.72*v:.2f}"/>')
    sx, sy = x + 12 + best[0] * cw + cw / 2 - 1, y + 20 + best[1] * ch + ch / 2
    star = " ".join(f"{sx + 9*math.cos(math.radians(a-90)):.1f},{sy + 9*math.sin(math.radians(a-90)):.1f} "
                    f"{sx + 3.9*math.cos(math.radians(a+36-90)):.1f},{sy + 3.9*math.sin(math.radians(a+36-90)):.1f}"
                    for a in range(0, 360, 72))
    s.append(f'<polygon points="{star}" fill="{RED}" stroke="#ffffff" stroke-width="1"/>')
    s.append(label(x + 104, y + 120, "period × duty", 10, 400, MUTED, "middle"))

    s.append(f'<rect x="700" y="542" width="682" height="30" rx="7" fill="{PAPER}" stroke="{LINE}"/>')
    s.append(label(714, 562, "ground truth for benchmarks · ABC-SMC · telegraph likelihood · treatment outcomes", 10, 400, INK))


CODE = """using Biomodelling, Random
# a resistance gene whose promoter switches more slowly than the cell divides
model = telegraph_model(k_on=0.005, k_off=0.005, k_tx=30.0, k_dm=1.0, k_tl=4.0, k_dp=0.2)
pert  = Perturbation(PiecewiseDose([0.0, 100.0], [0.0, 1.0]);
                     effects = [DeathHazard(h_max=0.5, EC50=0.3, protect=:protein, K=60.0)])
st    = PopulationSettings(dt=0.1, growth=ExponentialGrowth(log(2)/20), size_control=Sizer(2.0),
                           replication=Replication(0.5), control=FreeGrowth())
res   = simulate_population(model, x0, 1000, (0.0, 200.0); settings=st, perturbation=pert, rng=Xoshiro(1))

heritability(res, :protein; relation=:sisters)                          # expression memory across divisions
kill_curve(res); time_to_progression(res)                               # treatment outcomes
sequence(sample_cells(res; n=500).counts, SeqProtocol(capture=0.15)).Y   # synthetic scRNA-seq"""

CALLS = ("telegraph_model Perturbation PiecewiseDose DeathHazard PopulationSettings ExponentialGrowth Sizer "
         "Replication FreeGrowth simulate_population Xoshiro heritability kill_curve time_to_progression "
         "sequence sample_cells SeqProtocol log").split()


def highlight(line):
    """Colour a line of Julia: comment, keyword, call, number, everything else."""
    code, comment = line, ""
    hit = re.search(r"\s#\s", line)
    if hit:
        code, comment = line[:hit.start()], line[hit.start():]
    out, i = [], 0
    for m in re.finditer(r"[A-Za-z_][A-Za-z_0-9!]*|\d+\.?\d*", code):
        out.append(esc(code[i:m.start()]))
        tok = m.group(0)
        if tok == "using":
            out.append(f'<span style="color:{ORANGE}">{tok}</span>')
        elif tok in CALLS:
            out.append(f'<span style="color:{BLUE}">{tok}</span>')
        elif re.fullmatch(r"\d+\.?\d*", tok):
            out.append(f'<span style="color:{VIOLET}">{tok}</span>')
        else:
            out.append(esc(tok))
        i = m.end()
    out.append(esc(code[i:]))
    if comment:
        out.append(f'<span style="color:{MUTED};font-style:italic">{esc(comment)}</span>')
    return "".join(out)


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def add_panel_e(s):
    """e: a complete example. The code is an HTML overlay so the monospace metrics are exact."""
    s.append(panel(24, 596, 1392, 312, "e", "A complete experiment in twelve lines", AMBER))
    return s


def code_html():
    return ('<pre id="code">' + "\n".join(highlight(l) for l in CODE.split("\n")) + "</pre>")


HTML = """<!doctype html><html><head><meta charset="utf-8"><style>
@font-face {{ font-family: Inter; font-weight: 400; src: url(data:font/ttf;base64,{f400}) format('truetype'); }}
@font-face {{ font-family: Inter; font-weight: 600; src: url(data:font/ttf;base64,{f600}) format('truetype'); }}
@font-face {{ font-family: Inter; font-weight: 800; src: url(data:font/ttf;base64,{f800}) format('truetype'); }}
* {{ margin:0; padding:0; }}
body {{ width:{w}px; height:{h}px; background:#ffffff; overflow:hidden; position:relative; }}
#code {{ position:absolute; left:60px; top:638px; margin:0; font-size:15.6px; line-height:22.4px;
        font-family:"DejaVu Sans Mono","SFMono-Regular",Menlo,Consolas,monospace; color:#0b0b0b; white-space:pre; }}
svg {{ font-family: Inter, "Helvetica Neue", Arial, sans-serif; }}
text {{ dominant-baseline: alphabetic; }}
.mono {{ font-family: "DejaVu Sans Mono", "SFMono-Regular", Menlo, Consolas, monospace; }}
</style></head><body>
<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg">
<defs><marker id="arrowhead" markerWidth="7" markerHeight="7" refX="6" refY="3" orient="auto">
<path d="M0 0 L6.5 3 L0 6 z" fill="{ink2}"/></marker></defs>
{body}
</svg>
{code}
</body></html>"""


def make_html():
    w, h, s = make_figure1()
    add_panel_b(s)
    add_panel_c(s)
    add_panel_d(s)
    add_panel_e(s)
    return HTML.format(w=w, h=932, body="\n".join(s), ink2=INK2, code=code_html(),
                       f400=font(400), f600=font(600), f800=font(800))


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    out = os.path.join(FIG, "fig1_framework.html")
    with open(out, "w") as fh:
        fh.write(make_html())
    print("wrote", out)
