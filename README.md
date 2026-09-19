# Biomodelling.jl

[![CI](https://github.com/ayoublasri/Biomodelling.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ayoublasri/Biomodelling.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ayoublasri/Biomodelling.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ayoublasri/Biomodelling.jl)

Mechanistic stochastic simulation of gene regulatory networks inside **growing,
dividing and drug-treated cell populations**, with lineage tracking,
single-cell observation models and likelihood-free inference.

Authors: Ayoub Lasri and Marc Sturrock (Royal College of Surgeons in Ireland).
Version 1 of the framework was published in
[BMC Bioinformatics (2022)](https://doi.org/10.1186/s12859-022-04778-9); version
2.0 is a rewrite described in the accompanying preprint (see `paper/`).

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
  resistance proteins, growth inhibition, drug-induced rate changes, gene
  knockdown / overexpression in subsets of cells.
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

## Reproducing the paper

`paper/` contains the scripts that generate every figure of the v2.0 preprint
from fixed seeds (`julia --project=paper paper/run_all.jl`, then
`python paper/plot_all.py`).

## License

MIT. See `LICENSE`.
