# Changelog

## v2.0.0 (unreleased)

Breaking release. The package is reorganised around six layers: model, kernels,
cells and populations, perturbations, lineage, observation, plus inference and
generators. See `docs/plan/v2_implementation_and_paper_plan.md` for the design
and `docs/src/migration.md` for the migration guide.

### Added
- Every figure of the paper now sets in one typeface with one spacing scale, and
  Figure 6 gained the posterior predictive panel its legend already described.
- Supplementary Note 10 of the white paper: a claim-by-claim table separating what prior
  work established from what this work adds, with the comparison table in the main text
  now citing every tool it compares.
- `whitepaper/`: the version 2.0 white paper (PDF and DOCX) and its supplementary
  information, rebuilt from `paper/manuscript/` with the Nature citation style.
- Schedule optimisation: `AdaptiveDose` feedback schedules, finite-cycle `PulsedDose`,
  `daily_boluses`/`cycle_days` pharmacokinetic regimens, treatment outcomes
  (`time_to_progression`, `log_kill`, `net_growth_rate`, `extinction_probability`,
  `cumulative_dose`), `optimize_schedule`, state-dependent `GrowthInhibition`, `kin_pairs`.
- `paper/data/`: laboratory and clinical reference values (with provenance) used for calibration.
- `paper/`: manuscript sources, figure scripts (`paper/run_all.jl`, `paper/plot_all.py`), Python benchmark
  scripts and the built PDF/DOCX of the accompanying white paper.
- `ReactionModel` with `MassAction`, `Hill` and `Custom` kinetics, named
  parameters, volume scaling rules and promoter groups.
- Kernels `DirectSSA`, `TauLeap`, `HybridSSATau`, `AdaptiveTauLeap` sharing a
  preallocated `Workspace` and explicit random number generators.
- Cell and population layer: exponential growth with heterogeneity, `Sizer`,
  `Adder`, `AgeTimer` size control, binomial and beta-binomial partitioning,
  promoter inheritance, gene replication, `ConstantN`, `FreeGrowth` and
  `LogisticGrowth` population control, multithreaded population loop.
- Lineage recording (`LineageTable`) with heritability statistics,
  fluctuation tests and Newick export.
- Drug and perturbation layer: dose schedules, bolus pharmacokinetics,
  state-dependent death hazards, growth inhibition, fitness costs of cell
  states, rate modulation and gene perturbations; populations can start from
  per-cell snapshots.
- Observation models for scRNA-seq, smFISH and time-lapse reporters.
- Inference: ABC-SMC and the exact telegraph (Beta-Poisson) likelihood.
- Random gene regulatory network generator.
- CSV, Newick and AnnData-compatible HDF5 output (`HDF5` extension).
- Test suite with statistical validation against analytic distributions.

### Changed
- Minimum Julia version is 1.10.
- Continuous integration moved from Travis CI to GitHub Actions.

### Deprecated
- The v1 API (`Donne`, `ssa`, `tauleap`, `tauleapswitch`, `adaptive_tauleap`,
  `exponential_growth`, NamedTuple reactions) is kept as a thin compatibility
  layer and will be removed in v3.

### Removed
- `LsqFit` and `StringDistances` dependencies.
- The duplicated kernel files (`ssa2`, `ssa3`, `ssa_switch*`, `tauleapswitch2`,
  `non_negative_Poisson_tauleap`, `comp_*`, `compute_L`, `HO_reaction`).

## v0.3.0

Last release of the original framework (BMC Bioinformatics 2022).
