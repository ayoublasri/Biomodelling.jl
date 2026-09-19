# Reproducing the Biomodelling.jl 2.0 preprint

Every figure is generated from fixed seeds in two steps:

```bash
julia --project=paper -t 4 paper/run_all.jl        # simulations → paper/output/*.csv
python paper/plot_all.py                            # figures     → paper/figures/*.pdf|png
```

PIDC scores (Figure 5b) come from `paper/scripts/pidc_scores.jl`, which needs NetworkInference.jl in its own environment (`julia --project=paper/pidc paper/scripts/pidc_scores.jl`) because that package pins older dependencies. `paper/Project.toml` pins the Julia dependencies (add the package itself with
`Pkg.develop(path=".")` from the repository root); `paper/requirements.txt`
lists the Python plotting and benchmark tools.

Scripts (`paper/scripts/`):

| Script | Figure | Content |
|---|---|---|
| `fig2_engine.jl` | 2 | exact distributions, kernel agreement, JumpProcesses cross-check, runtime |
| `fig3_growth.jl` | 3 | lineage traces, size scaling, heritability, memory timescales, lineage vs population noise, copy number |
| `fig4_persisters.jl` | 4 | kill curves, decay rate vs dose, fate correlations, clone diversity, schedules, memory disruption, MGMT |
| `fig5_benchmarks.jl` + `benchmarks_python.py` | 5 | memory genes, GRN inference, imputation regression, perturbation baselines |
| `fig6_inference.jl` | 6 | ABC vs exact likelihood, division bias, drug parameter recovery |

The manuscript sources are in `paper/manuscript/`.
