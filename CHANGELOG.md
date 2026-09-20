# Changelog

## v2.0.0 (unreleased)

Breaking release. The package is reorganised around six layers: model, kernels,
cells and populations, perturbations, lineage, observation, plus inference and
generators. See `docs/plan/v2_implementation_and_paper_plan.md` for the design
and `docs/src/migration.md` for the migration guide.

### Added
- `MotherMachine` population control: at every division one daughter is kept at
  random and the other discarded, so `N0` founders give `N0` independent
  single-lineage traces. With it, the single-lineage and population-snapshot
  distributions can be generated from the same model and compared.
- `AgeTimer` accepts gamma-distributed interdivision times (`dist = :gamma`),
  covering the Erlang family of the exactly solvable population models and, at
  `cv = 1`, the memoryless timer.
- `simulate_population` accepts `age0`, the founders' ages.
- Validation of the population and lineage layers against the exact stationary
  solutions that exist for this model class (Beentjes, Perez-Carrasco & Grima
  2020; Jia & Grima 2023), in both lineage and snapshot modes, for constitutive,
  bursty, telegraph and volume-scaled-with-replication kinetics
  (`paper/scripts/fig9_validation.jl` and `paper/scripts/exact_solutions.py`,
  reported as Supplementary Note 14). Single lineages match the exact law of
  their mode to Kolmogorov-Smirnov distances of 0.0023 to 0.0054 against 99%
  critical values of 0.0081 to 0.0115; population snapshots, whose cells are not
  independent, agree in the mean to within 2.5 standard errors; and the same
  samples scored against the law of the other mode give distances at least six
  times larger. The closed-form lineage and population means and variances are
  checked on every run of the test suite.
- `CycleSensitivity`: cell-cycle dependence of the drug death hazard, for agents
  whose lesions are converted into death during replication. The hazard is scaled by
  `baseline + (1 - baseline)·exp(-((φ - center)/width)²/2)` at cycle progress `φ`.
- `SuicideConsumption`: stoichiometric consumption of a protective protein by the
  drug (one molecule per lesion repaired, at rate `k·d·V·c/(c + K_m)`), so that
  depletion follows the cumulative exposure rather than the peak concentration.
  This replaces the dose-scaled first-order degradation used for MGMT in Fig. 8.
- Cell-cycle robustness runs for the persister and melanoma case studies
  (`fig4_persisters.jl cycle`, `fig8_schedules.jl melanoma_cycle`), reported as
  Supplementary Note 13 and Supplementary Figure 5.
- Identifiability and robustness analysis of the calibration (`fig7_calibration.jl profile`):
  conditional parameter profiles, a memory/resistant-fraction slice and an integration-step
  check, reported as Supplementary Note 11 and Supplementary Figure 3.
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

### Security
- Upgraded the `video/` build tooling from Remotion 4.0.416 to 4.0.526, clearing
  8 advisories (7 high) in transitive build dependencies: `extract-zip` symlink
  path traversal and arbitrary file write, `webpack` buildHttp SSRF via allowlist
  bypass, and `ws` memory disclosure and DoS. `npm audit` now reports none.
- Every GitHub Actions workflow declares least-privilege `permissions:`.
  CompatHelper no longer runs hourly on Julia 1.2 (a version that cannot load this
  package) via an unpinned `@latest` action; it runs daily on Julia 1 with a pinned
  action. TagBot is triggered by the registrator's comment rather than polling
  hourly, and is gated on the actor.
- Added `.github/dependabot.yml` (weekly npm and Actions updates), `SECURITY.md`
  with a private reporting route, and a CI job that fails on any moderate or worse
  npm advisory.
- Stopped committing the Inter font binaries; the video build fetches them, matching
  how the paper figures already worked.

### Changed
- Minimum Julia version is 1.10.
- Manuscript revised in response to peer review: the memory-disruption
  prediction is reconciled with the two published pretreatment experiments that
  bear on it, the temozolomide consumption claim is softened against the
  measured depletion time-courses, an explicit limitation is added for the kin
  correlations that independent cycle noise cannot generate, the capability
  table and the novelty claim are corrected, and every reported error metric is
  given one definition. The cell-cycle gate comparison gains an arm at matched
  cycle-averaged mean hazard, the schedule scans gain replicate seeds, and the
  melanoma memory assumption gains a sensitivity scan.
- Continuous integration moved from Travis CI to GitHub Actions.
- Manuscript revised after a multi-dimension review (citations, numbers, biology,
  figures, methods-vs-code, prose): corrected the Fig. 3c heritability description,
  the memory-disruption and imputation results, the hybrid-kernel exactness wording
  and the Fig. 6d and Fig. 7b axes; cited the primary Das Thakur et al. (2013, Nature)
  xenograft study; labelled the MGMT fractions and consumption model, the RECIST
  analogue and the temozolomide plasma half-life as the approximations they are;
  aligned the parameter table (panels 8a-e / 8f-h, 150 particles, 500 founders,
  dt and profile settings) with the scripts; and fixed the Fig. 1e code snippet
  (`initial_state`, `log_kill`).

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
