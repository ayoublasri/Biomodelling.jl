# Biomodelling.jl

Mechanistic stochastic simulation of gene regulatory networks inside growing,
dividing and drug-treated cell populations.

Biomodelling.jl couples exact and approximate stochastic simulation of reaction
networks (Gillespie's direct method, tau-leaping, hybrid and adaptive schemes)
with the physiology that shapes single-cell data: exponential growth,
volume-scaled transcription, gene replication, division with binomial
partitioning of molecules and inheritance of promoter states, population
control, and drug-induced death. It records the full division tree, so
heritable ("memory") expression states, lineage correlations and
Luria-Delbrück fluctuation statistics can be measured directly; it converts
true molecule counts into scRNA-seq, smFISH or time-lapse observations; and it
provides likelihood-free inference (ABC-SMC) and the exact telegraph-model
likelihood.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/ayoublasri/Biomodelling.jl")
```

Julia 1.10 or later is required. Optional: `HDF5` for AnnData (`.h5ad`) export.

## Ten-line example

```julia
using Biomodelling, Random

model = telegraph_model(k_on = 0.01, k_off = 0.01, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
x0    = initial_state(model; G_off = 1)
st    = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20),
                           size_control = Sizer(2.0; cv = 0.05), replication = Replication(0.5))
res   = simulate_population(model, x0, 500, (0.0, 200.0); settings = st, rng = Xoshiro(1))

heritability(res, :protein; relation = :mother_daughter)   # memory across divisions
sn = sample_cells(res; n = 300)                            # a snapshot of 300 cells
Y  = sequence(sn.counts, SeqProtocol(capture = 0.15)).Y    # synthetic scRNA-seq counts
```

## Layers

| Layer | Main types and functions |
|---|---|
| Model | [`ReactionModel`](@ref), [`Reaction`](@ref), [`MassAction`](@ref), [`Hill`](@ref), [`Custom`](@ref), builders |
| Kernels | [`DirectSSA`](@ref), [`TauLeap`](@ref), [`HybridSSATau`](@ref), [`AdaptiveTauLeap`](@ref), [`simulate`](@ref), [`ensemble_final`](@ref) |
| Cells and populations | [`PopulationSettings`](@ref), [`simulate_population`](@ref), [`ExponentialGrowth`](@ref), [`Sizer`](@ref), [`Adder`](@ref), [`AgeTimer`](@ref), [`BinomialPartition`](@ref), [`Replication`](@ref), [`ConstantN`](@ref), [`FreeGrowth`](@ref), [`LogisticGrowth`](@ref) |
| Perturbations | [`Perturbation`](@ref), dose schedules, [`DeathHazard`](@ref), [`GrowthInhibition`](@ref), [`GrowthCost`](@ref), [`CycleSensitivity`](@ref), [`SuicideConsumption`](@ref), [`RateModulation`](@ref), [`GenePerturbation`](@ref) |
| Lineage | [`LineageTable`](@ref), [`heritability`](@ref), [`lineage_autocorrelation`](@ref), [`memory_timescale`](@ref), [`fluctuation_test`](@ref), [`clonal_variance_scores`](@ref), [`newick`](@ref) |
| Observation | [`SeqProtocol`](@ref), [`sequence`](@ref), [`smfish`](@ref), [`timelapse`](@ref), [`sample_cells`](@ref) |
| Inference | [`abc_smc`](@ref), [`telegraph_pmf`](@ref), [`fit_telegraph`](@ref), [`moment_summaries`](@ref) |
| Generators and IO | [`random_grn`](@ref), CSV writers, [`write_h5ad`](@ref) |

## Citing

If you use Biomodelling.jl, please cite the software itself using the metadata in
`CITATION.cff` at the repository root, together with the version 2.0 white paper,
*Mechanistic simulation of heritable expression states, cell division and drug
response in single-cell populations* (in the `whitepaper/` directory of the
repository), and the version 1 article (Lasri, Shahrezaei and Sturrock, *BMC
Bioinformatics* 23:236, 2022).
