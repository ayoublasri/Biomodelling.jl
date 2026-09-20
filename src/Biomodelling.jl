"""
    Biomodelling

Mechanistic stochastic simulation of gene regulatory networks inside growing,
dividing and drug-treated cell populations, with lineage tracking, observation
models for single-cell measurements, and likelihood-free inference.

The package is organised in layers:

* `model`   : `ReactionModel`, `Reaction`, kinetics (`MassAction`, `Hill`, `Custom`)
* `kernels` : `DirectSSA`, `TauLeap`, `HybridSSATau`, `AdaptiveTauLeap`, `simulate`
* `cells`   : growth, size control, division, partitioning, `simulate_population`
* `perturb` : dose schedules and drug effects
* `lineage` : `LineageTable`, heritability statistics, fluctuation tests
* `observe` : scRNA-seq, smFISH and time-lapse observation models
* `infer`   : summary statistics, ABC-SMC, telegraph likelihood
"""
module Biomodelling

using Random
using Random: Xoshiro, randexp
using Statistics
using LinearAlgebra
using Printf
using Distributions
using StatsBase: StatsBase, sample, countmap, cor, autocor, Weights
using PoissonRandom: pois_rand
using SpecialFunctions: loggamma
using QuadGK: quadgk
using Optim
using DelimitedFiles

# ---- model -----------------------------------------------------------------
include("model/reactions.jl")
include("model/model.jl")
include("model/builders.jl")
# ---- kernels ---------------------------------------------------------------
include("kernels/workspace.jl")
include("kernels/propensities.jl")
include("kernels/direct_ssa.jl")
include("kernels/tauleap.jl")
include("kernels/adaptive_tauleap.jl")
include("kernels/simulate.jl")
# ---- cells -----------------------------------------------------------------
include("cells/cell.jl")
include("cells/growth.jl")
include("cells/division.jl")
include("cells/partitioning.jl")
# ---- perturbations ---------------------------------------------------------
include("perturb/schedules.jl")
include("perturb/effects.jl")
# ---- lineage table (needed by the population loop) -------------------------
include("lineage/tree.jl")
# ---- population loop -------------------------------------------------------
include("cells/population.jl")
include("perturb/optimize.jl")
# ---- lineage statistics ----------------------------------------------------
include("lineage/heritability.jl")
include("lineage/fluctuation_test.jl")
# ---- observation -----------------------------------------------------------
include("observe/scrnaseq.jl")
include("observe/timelapse.jl")
# ---- inference -------------------------------------------------------------
include("infer/summaries.jl")
include("infer/abc_smc.jl")
include("infer/telegraph_likelihood.jl")
# ---- generators and io -----------------------------------------------------
include("generators/random_grn.jl")
include("io/tables.jl")
include("io/newick.jl")
# ---- v1 compatibility ------------------------------------------------------
include("compat_v1.jl")

# model
export Species, Reaction, ReactionModel, MassAction, Hill, Custom
export nspecies, nreactions, nparams, speciesnames, paramnames, speciesindex, paramindex
export set_params, initial_state, propensities, stoichiometry
export telegraph_model, birth_death_model, two_stage_model, bursty_protein_model
# kernels
export AbstractKernel, DirectSSA, TauLeap, HybridSSATau, AdaptiveTauLeap, Workspace
export simulate, ensemble_final, Trajectory
# cells
export Cell, ExponentialGrowth, Sizer, Adder, AgeTimer, BinomialPartition, BetaBinomialPartition
export Replication, ConstantN, FreeGrowth, LogisticGrowth, PopulationSettings, PopulationResult
export simulate_population, snapshot, final_snapshot, popsize, concentrations
# perturbation
export DoseSchedule, ConstantDose, PulsedDose, PiecewiseDose, BolusPK, FunctionDose, AdaptiveDose, dose, daily_boluses, cycle_days, cumulative_dose
export net_growth_rate, log_kill, time_to_progression, extinction_probability, optimize_schedule, ScheduleOptimum
export DrugEffect, DeathHazard, GrowthInhibition, GrowthCost, RateModulation, GenePerturbation, Perturbation
export CycleSensitivity, SuicideConsumption
# lineage
export LineageTable, children, sister_pairs, cousin_pairs, kin_pairs, mother_daughter_pairs, lineage_of, follow_lineage
export heritability, lineage_autocorrelation, memory_timescale, noise_decomposition
export fluctuation_test, clonal_variance_scores, covariance_eigenspectrum, powerlaw_tail_exponent
# observation
export SeqProtocol, sequence, smfish, timelapse, sample_cells
# inference
export moment_summaries, summary_distance, abc_smc, ABCResult, posterior_mean, posterior_median, posterior_quantile, credible_interval
export telegraph_pmf, telegraph_loglik, fit_telegraph
# generators & io
export random_grn, write_counts_csv, write_lineage_csv, write_metadata_csv, newick, write_h5ad
# v1 compatibility
export Donne, ssa, tauleap, tauleapswitch, adaptive_tauleap, exponential_growth, growth_estimate

end # module
