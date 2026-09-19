# Reproducing the Biomodelling.jl 2.0 preprint

Every figure is generated from fixed seeds in two steps:

```bash
julia --project=paper -t 4 paper/run_all.jl        # simulations → paper/output/*.csv
python paper/plot_all.py                            # figures     → paper/figures/*.pdf|png
```

`paper/Project.toml` pins the Julia dependencies (add the package itself with
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
| `fig7_calibration.jl profile` | S3 | conditional parameter profiles, the memory/resistant-fraction slice and the integration-step check |
| `fig7_calibration.jl` | 7 | calibration of the persister model to the U2OS cisplatin fate data of Iyer et al. (2025) by simulated moments; validation on the held-out concentration, kin correlations, death-time invariance and kill-curve shape |
| `fig8_schedules.jl` | 8 | continuous, intermittent (SWOG S1320) and adaptive schedules for three resistance mechanisms, schedule optimiser maps, temozolomide regimens of RTOG 0525 with MGMT stable or consumed, fractionation at equal cumulative dose |

The manuscript sources are in `paper/manuscript/`.
