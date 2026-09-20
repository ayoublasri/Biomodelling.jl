# Biomodelling.jl

[![CI](https://github.com/ayoublasri/Biomodelling.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ayoublasri/Biomodelling.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ayoublasri/Biomodelling.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ayoublasri/Biomodelling.jl)

Mechanistic stochastic simulation of gene regulatory networks inside **growing,
dividing and drug-treated cell populations**, with lineage tracking,
single-cell observation models and likelihood-free inference.

Author: Ayoub Lasri. Version 2.0 was
authored with [Claude](https://claude.ai) (Anthropic) under the author's
direction. Version 1 of the framework was developed with Marc Sturrock and
Vahid Shahrezaei and published in
[BMC Bioinformatics (2022)](https://doi.org/10.1186/s12859-022-04778-9).

**White paper:** [Mechanistic simulation of heritable expression states, cell
division and drug response in single-cell populations](whitepaper/Biomodelling-jl-2.0-white-paper.pdf)
([supplementary information](whitepaper/Biomodelling-jl-2.0-supplementary-information.pdf)).

## What it does

* **Reaction networks** with mass-action, Hill (activation, inhibition,
  combinatorial) and custom kinetics, named parameters, and volume scaling
  rules that keep concentrations consistent as cells grow.
* **Stochastic kernels**: Gillespie's direct method, fixed-step tau-leaping, a
  hybrid tau-leap/SSA scheme and adaptive tau-leaping (Cao, Gillespie and
  Petzold), all reproducible from explicit random number generators.
* **Cells and populations**: exponential growth with heterogeneity, sizer /
  adder / timer division control, binomial or beta-binomial partitioning,
  promoter-state inheritance, gene replication, constant-size, free-growth or
  logistic population control, multithreaded and deterministic.
* **Drug and perturbation layer**: dose schedules (constant, pulsed, piecewise,
  one-compartment pharmacokinetics), state-dependent death hazards protected by
  resistance proteins, growth inhibition, fitness costs of resistant states,
  drug-induced rate changes, gene knockdown / overexpression in subsets of cells.
* **Lineage**: the complete division tree; mother-daughter, sister and cousin
  correlations; lineage autocorrelation and memory timescales; lineage versus
  population noise; Luria-Delbrück fluctuation tests; MemorySeq-style clonal
  scores; Newick export.
* **Observation models**: scRNA-seq (capture efficiency, depth, dropout,
  batches), smFISH and time-lapse reporters; CSV and AnnData (`.h5ad`) output.
* **Inference**: ABC-SMC over any simulation, and the exact Beta-Poisson
  likelihood of the telegraph model.
* **Generators**: random regulatory networks with known ground truth.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/ayoublasri/Biomodelling.jl")
```

Requires Julia 1.10 or later. `using HDF5` enables `write_h5ad`.

## Quick start

```julia
using Biomodelling, Random

# a resistance gene with slow promoter switching (heritable expression states)
model = telegraph_model(k_on = 0.005, k_off = 0.005, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
x0    = initial_state(model; G_off = 1)

# growing, dividing population; drug applied from t = 100 kills low-expressing cells
pert  = Perturbation(PiecewiseDose([0.0, 100.0], [0.0, 1.0]);
                     effects = [DeathHazard(h_max = 0.5, EC50 = 0.3, protect = :protein, K = 60.0, q = 4.0)])
st    = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20), size_control = Sizer(2.0; cv = 0.05),
                           replication = Replication(0.5), control = FreeGrowth(max_cells = 20_000))
res   = simulate_population(model, x0, 1000, (0.0, 200.0); settings = st, perturbation = pert, rng = Xoshiro(1))

res.popsize                                                # kill curve and regrowth
heritability(res, :protein; relation = :mother_daughter)   # expression memory across divisions
sn = sample_cells(res; n = 500)                            # snapshot of survivors
Y  = sequence(sn.counts, SeqProtocol(capture = 0.15)).Y    # synthetic scRNA-seq counts
```

See the documentation (`docs/`) for tutorials on single cells, populations,
drug treatment, synthetic data and inference, and `docs/src/migration.md` for
the mapping from the v1 API (which still works through a deprecated
compatibility layer).

## White paper

`whitepaper/` holds the v2.0 white paper,
[*Mechanistic simulation of heritable expression states, cell division and drug
response in single-cell populations*](whitepaper/Biomodelling-jl-2.0-white-paper.pdf),
with its [supplementary information](whitepaper/Biomodelling-jl-2.0-supplementary-information.pdf)
and an [editable version](whitepaper/Biomodelling-jl-2.0-white-paper.docx). It
reports, with all code and seeds in this repository:

* validation of the kernels against exact stationary laws and against
  JumpProcesses.jl;
* cell-size scaling, partitioning noise, lineage-versus-population statistics
  and cell-cycle-dependent bursting arising from the physiology alone;
* a single resistance gene with slow promoter switching reproducing the
  signatures of drug-tolerant persisters (dose-dependent decay with unchanged
  single-cell timing, correlated fates of related cells, barcode diversity
  preserved or collapsed, loss of resistance when memory is disrupted during
  exposure);
* ground-truth benchmarks for memory-gene detection, network inference,
  imputation and perturbation prediction;
* calibration to published time-lapse measurements of cisplatin-treated cells
  and validation on a held-out concentration, on lineage correlations and on
  the dose-invariance of single-cell timing;
* schedule optimisation reproducing the opposite outcomes of intermittent
  dosing in the SWOG S1320 melanoma trial and in patient-derived xenografts,
  and evaluating standard against dose-dense temozolomide (RTOG 0525).

The paper sources, figure scripts and reference data are in `paper/`; every
figure is regenerated from fixed seeds with `julia --project=paper
paper/run_all.jl` followed by `python paper/plot_all.py`, and the document is
rebuilt with `bash paper/manuscript/build.sh`.

## Announcement assets

`whitepaper/social/` holds the announcement card (`linkedin-card.png`, 1200 x
1200) with the script that renders it from the schedule results, and the post
text.

## Citation

If you use Biomodelling.jl, please cite the white paper above and the version 1
article: Lasri, A., Shahrezaei, V. and Sturrock, M. Benchmarking imputation
methods for network inference using a novel method of synthetic scRNA-seq data
generation. *BMC Bioinformatics* **23**, 236 (2022).
<https://doi.org/10.1186/s12859-022-04778-9>

## License

MIT. See `LICENSE`.

## Dose and schedule optimisation

The perturbation layer includes clinical-style regimens (`PulsedDose` with a
finite number of cycles, `daily_boluses` with one-compartment pharmacokinetics,
`cycle_days`), feedback schedules (`AdaptiveDose`, adaptive therapy), treatment
outcomes (`time_to_progression`, `log_kill`, `net_growth_rate`,
`extinction_probability`, `cumulative_dose`) and a bounded optimiser with common
random numbers (`optimize_schedule`) that minimises any outcome over schedule
parameters. `paper/scripts/fig7_calibration.jl` calibrates the persister model
to published time-lapse data and validates it on held-out drug concentrations;
`paper/scripts/fig8_schedules.jl` compares continuous, intermittent (SWOG S1320)
and adaptive dosing and the standard against dose-dense temozolomide regimens
(RTOG 0525) and searches for better schedules. See the documentation tutorial
"Dose and schedule optimisation".
