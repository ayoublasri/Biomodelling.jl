"""Render the LinkedIn card for Biomodelling.jl 2.0.

    python build_card.py                     # writes card.html next to this file
    chromium --headless --disable-gpu --hide-scrollbars --force-device-scale-factor=2 \
        --window-size=1200,1340 --virtual-time-budget=4000 --screenshot=raw.png card.html
    # then crop the top 2400x2400 pixels of raw.png

Numbers are the mean weeks to loss of control from paper/output/fig8a_melanoma_schedules.csv.
"""
import base64, os, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
FONTS = {  # Inter, served by Google Fonts; downloaded once into fonts/
    "400": "https://fonts.gstatic.com/s/inter/v20/UcCO3FwrK3iLTeHuS_nVMrMxCp50SjIw2boKoduKmMEVuLyfMZg.ttf",
    "600": "https://fonts.gstatic.com/s/inter/v20/UcCO3FwrK3iLTeHuS_nVMrMxCp50SjIw2boKoduKmMEVuGKYMZg.ttf",
    "800": "https://fonts.gstatic.com/s/inter/v20/UcCO3FwrK3iLTeHuS_nVMrMxCp50SjIw2boKoduKmMEVuDyYMZg.ttf",
}
def font(w):
    """Base64 of Inter at weight w, downloading it on first use; '' falls back to the system sans."""
    path = os.path.join(HERE, "fonts", f"Inter-{w}.ttf")
    if not os.path.exists(path):
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            urllib.request.urlretrieve(FONTS[str(w)], path)
        except Exception as e:
            print(f"  [Inter {w} unavailable ({e}); falling back to the system sans-serif]")
            return ""
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

SER = [("Continuous dosing", "#2a78d6"), ("Intermittent (trial schedule)", "#eb6834"), ("Adaptive therapy", "#1baf7a")]
GROUPS = [("Resistance is free", "no fitness cost", [23.4, 26.1, 26.7], [0,0,0]),
          ("Resistance costs growth", "50% slower without drug", [52.0, 52.0, 52.0], [1,1,1]),
          ("Resistance is partial", "resistant cells slowed by drug", [52.0, 48.6, 42.7], [1,0,0])]

W, H = 1040, 400          # plot area
L, R, T, B = 8, 8, 40, 64  # margins inside the svg
pw, ph = W - L - R, H - T - B
YMAX = 60.0
def y(v): return T + ph * (1 - v / YMAX)

gap, bar_gap = 64, 6
gw = (pw - gap * (len(GROUPS) - 1)) / len(GROUPS)
bw = (gw - bar_gap * (len(SER) - 1)) / len(SER)

bars, labels, ticks = [], [], []
for gi, (name, sub, vals, cens) in enumerate(GROUPS):
    gx = L + gi * (gw + gap)
    for si, (v, c) in enumerate(zip(vals, cens)):
        x = gx + si * (bw + bar_gap)
        top, base = y(v), y(0)
        bars.append(f'<path d="M{x:.1f} {base:.1f} L{x:.1f} {top+4:.1f} Q{x:.1f} {top:.1f} {x+4:.1f} {top:.1f} '
                    f'L{x+bw-4:.1f} {top:.1f} Q{x+bw:.1f} {top:.1f} {x+bw:.1f} {top+4:.1f} L{x+bw:.1f} {base:.1f} Z" fill="{SER[si][1]}"/>')
        txt = f"{v:.0f}+" if c else f"{v:.0f}"
        labels.append(f'<text x="{x+bw/2:.1f}" y="{top-12:.1f}" class="val">{txt}</text>')
    labels.append(f'<text x="{gx+gw/2:.1f}" y="{y(0)+30:.1f}" class="glab">{name}</text>')
    labels.append(f'<text x="{gx+gw/2:.1f}" y="{y(0)+50:.1f}" class="gsub">{sub}</text>')
for v in (0, 20, 40, 60):
    ticks.append(f'<line x1="{L}" y1="{y(v):.1f}" x2="{L+pw}" y2="{y(v):.1f}" stroke="#E6E5E1" stroke-width="1"/>')
cens_y = y(52)
ticks.append(f'<line x1="{L}" y1="{cens_y:.1f}" x2="{L+pw}" y2="{cens_y:.1f}" stroke="#8A8985" stroke-width="1.5" stroke-dasharray="6 5"/>')
ticks.append(f'<text x="{L}" y="{cens_y-11:.1f}" class="cens">end of 52-week follow-up</text>')

legend = "".join(
    f'<span class="lg"><i style="background:{c}"></i>{n}</span>' for n, c in SER)

html = f"""<!doctype html><html><head><meta charset="utf-8"><style>
@font-face {{ font-family: Inter; font-weight: 400; src: url(data:font/ttf;base64,{font(400)}) format('truetype'); }}
@font-face {{ font-family: Inter; font-weight: 600; src: url(data:font/ttf;base64,{font(600)}) format('truetype'); }}
@font-face {{ font-family: Inter; font-weight: 800; src: url(data:font/ttf;base64,{font(800)}) format('truetype'); }}
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ width:1200px; height:1200px; overflow:hidden; background:#FCFCFB; font-family:Inter, sans-serif; color:#0B0B0B;
       -webkit-font-smoothing:antialiased; display:flex; flex-direction:column; }}
.rule {{ height:10px; background:linear-gradient(90deg,#2a78d6 0 33.3%,#eb6834 33.3% 66.6%,#1baf7a 66.6% 100%); }}
.wrap {{ padding:52px 80px 0; flex:1; display:flex; flex-direction:column; }}
.eyebrow {{ font-size:19px; font-weight:600; letter-spacing:.14em; text-transform:uppercase; color:#52514E; }}
h1 {{ font-size:57px; line-height:1.1; font-weight:800; letter-spacing:-.022em; margin-top:16px; }}
h1 em {{ font-style:normal; color:#2a78d6; }}
.sub {{ font-size:23px; line-height:1.45; color:#52514E; margin-top:20px; max-width:1000px; font-weight:400; }}
.chart-title {{ font-size:20px; font-weight:600; color:#0B0B0B; margin-top:34px; }}
.note {{ font-size:17px; color:#8A8985; font-weight:400; }}
.legend {{ display:flex; gap:28px; margin-top:12px; }}
.lg {{ display:flex; align-items:center; gap:9px; font-size:18px; color:#52514E; font-weight:500; }}
.lg i {{ width:15px; height:15px; border-radius:4px; display:inline-block; }}
svg {{ margin-top:14px; }}
.val {{ font-size:22px; font-weight:800; fill:#0B0B0B; text-anchor:middle; }}
.glab {{ font-size:21px; font-weight:600; fill:#0B0B0B; text-anchor:middle; }}
.gsub {{ font-size:18px; font-weight:400; fill:#8A8985; text-anchor:middle; }}
.cens {{ font-size:16px; fill:#8A8985; text-anchor:start; }}
.ylab {{ font-size:18px; color:#8A8985; }}
.proof {{ display:flex; gap:22px; margin-top:auto; margin-bottom:30px; }}
.pc {{ flex:1; background:#F4F3F0; border-radius:14px; padding:20px 22px; }}
.pc b {{ display:block; font-size:19px; font-weight:600; color:#0B0B0B; line-height:1.3; }}
.pc span {{ display:block; font-size:17px; color:#52514E; margin-top:7px; line-height:1.35; }}
.foot {{ margin-top:auto; height:126px; flex:0 0 126px; display:flex; align-items:center; justify-content:space-between;
        border-top:1px solid #E6E5E1; margin-left:-80px; margin-right:-80px; padding:0 80px; }}
.brand {{ font-size:28px; font-weight:800; letter-spacing:-.02em; }}
.brand span {{ color:#52514E; font-weight:600; }}
.tag {{ font-size:19px; color:#52514E; margin-top:8px; }}
.url {{ font-size:19px; color:#2a78d6; font-weight:600; text-align:right; line-height:1.6; }}
</style></head><body>
<div class="rule"></div>
<div class="wrap">
  <div class="eyebrow">Open source · Julia · in silico single cells</div>
  <h1>The best drug schedule flips<br><em>on one property of the cells</em></h1>
  <div class="sub">A simulated melanoma carrying rare, inherited drug-tolerant cells. One model reproduces both a
  clinical trial and the xenograft result that contradicted it, once you set how fast resistant cells grow without the drug.</div>
  <div class="chart-title">Weeks until the tumour regrows past its starting size <span class="note">· + = still controlled when follow-up ended</span></div>
  <div class="legend">{legend}</div>
  <svg width="{W}" height="{H}" viewBox="0 0 {W} {H}">
    {''.join(ticks)}
    {''.join(bars)}
    {''.join(labels)}
  </svg>
  <div class="proof">
    <div class="pc"><b>Calibrated on real data</b><span>Fates of cells tracked by time-lapse under three cisplatin doses</span></div>
    <div class="pc"><b>Validated out of sample</b><span>Held-out dose, lineage correlations, dose-invariant death times</span></div>
    <div class="pc"><b>Reproducible by anyone</b><span>Fixed seeds, 203 tests, white paper and code under MIT</span></div>
  </div>
  <div class="foot">
    <div>
      <div class="brand">Biomodelling.jl <span>2.0</span></div>
      <div class="tag">Cells that grow, divide, inherit expression states and respond to drugs</div>
    </div>
    <div class="url">github.com/ayoublasri/Biomodelling.jl<br>White paper + code, MIT licence</div>
  </div>
</div>
</body></html>"""
open(os.path.join(HERE, "card.html"), "w").write(html)
print("wrote card.html", len(html))
