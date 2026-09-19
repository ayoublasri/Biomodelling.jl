# Biomodelling.jl 2.0: implementation plan and bioRxiv paper plan

Status: DRAFT FOR REVIEW. Nothing in the package has been changed yet. Every section below is a proposal; items marked **DECISION** need your answer before the corresponding work starts.

Prepared on the branch `claude/single-cell-drug-discovery-opportunities-7va73w`.

---

## 0. What you are approving

1. Reposition Biomodelling.jl from "synthetic scRNA-seq generator for imputation benchmarking" (the 2022 paper) to **the mechanistic simulator of heritable gene-expression states in growing, dividing, drug-treated cell populations**, aimed at non-genetic drug tolerance (persisters), lineage-resolved single-cell data, and benchmark generation.
2. Rebuild the package as v2.0 with a breaking but backward-shimmed API: one model type, one consolidated set of stochastic kernels, and five new layers (cells and division, drug and perturbation, lineage, observation, inference).
3. Write a bioRxiv research article (not an application note) with six main figures whose central validation is against (a) exact analytical results, (b) an independent simulator, (c) the authors' own 2020 and 2022 published results, and (d) quantitative findings from four recent persister and memory studies.
4. Deliver a reproducibility package (`paper/` scripts, seeds, pinned environment) so every figure regenerates from one command.

---

## 1. Decisions I need from you

| # | Decision | My recommendation | Why |
|---|----------|-------------------|-----|
| D1 | Package name and version | Keep `Biomodelling.jl`, release as **v2.0.0** (semver-breaking) | Registered name, continuity with the 2022 paper and its citations; a rename throws that away. |
| D2 | Simulation engine strategy | **Native kernels, consolidated and rewritten, plus a Catalyst.jl import bridge**; not a full SciML rebuild | Population simulations run thousands of independent cells for short intervals between divisions; a preallocated native kernel is faster and simpler than re-creating JumpProcesses problems per cell per step. The bridge (`ReactionSystem -> Biomodelling model`) gives interoperability; the cross-check against JumpProcesses.jl gives credibility. |
| D3 | Flagship biology for the drug layer | A **generic "resistance gene" persister model** calibrated to the published quantitative results of Iyer, Granada and Chakrabarti (cisplatin, PLoS Comput Biol 2025), plus schedule optimisation compared to Corigliano et al. (PRX Life 2025), plus a **re-implementation of your 2020 MGMT/temozolomide model** as a continuity check | Uses numbers already in the literature (no data access needed), and ties directly to your own prior work. If you have RCSI single-cell or time-lapse data under drug, say so: it would become the primary calibration target. |
| D4 | Inference scope | **ABC-SMC implemented in Julia** (likelihood-free, summary statistics), with exact likelihoods for the telegraph class used as a reference. Neural posterior estimation (Python `sbi`) only as a documented export path, not a deliverable | Keeps the package Julia-only and testable; NPE would add a heavy Python dependency and training cost. |
| D5 | Python dependencies for benchmarks | Allow Python tools **only in `paper/` benchmark scripts** (arboreto for GENIE3/GRNBoost2, MAGIC, kNN-smoothing); the package itself stays pure Julia (PIDC via NetworkInference.jl) | Same imputation and GRN methods as the 2022 paper where practical; package stays lightweight. |
| D6 | Paper venue and format | bioRxiv, category **Systems Biology**, research article (~6,000 words, 6 main figures, supplement); later journal target PLoS Computational Biology or NAR Genomics and Bioinformatics | Application-note length cannot carry the validation story. |
| D7 | Authors and affiliations | Placeholders: A. Lasri, M. Sturrock, with V. Shahrezaei as optional co-author (2022 co-author) | Your call. |
| D8 | License | Add an **MIT** license file (the repository currently has no LICENSE file) | Required for bioRxiv code availability and Julia registry hygiene. |
| D9 | Julia baseline | **Julia 1.10 (LTS)** minimum; CI on 1.10 and latest stable | Current `Project.toml` says 1.3 and CI is on defunct Travis. |
| D10 | Extra optional features | Nascent/mature RNA split (RNA velocity ground truth) and cell death by apoptosis-independent aging: **include nascent/mature as an option, skip the rest** | Cheap to add, widens benchmark use; the rest is scope creep. |

---

## 2. Package changes for v2.0

### 2.1 Repository restructure

```
Biomodelling.jl/
  Project.toml            # julia 1.10, compat bounds, version 2.0.0
  LICENSE                 # MIT (new)
  README.md               # rewritten quick start
  CHANGELOG.md            # new
  .github/workflows/CI.yml, Docs.yml, CompatHelper.yml, TagBot.yml
  src/
    Biomodelling.jl       # module, exports
    model/    reactions.jl  regulation.jl  model.jl  catalyst_bridge.jl  compat_v1.jl
    kernels/  direct_ssa.jl  tauleap.jl  adaptive_tauleap.jl  hybrid.jl  propensities.jl
    cells/    cell.jl  growth.jl  division.jl  partitioning.jl  population.jl
    perturb/  schedules.jl  pharmacokinetics.jl  effects.jl  gene_perturbations.jl
    lineage/  tree.jl  barcodes.jl  heritability.jl  fluctuation_test.jl
    observe/  scrnaseq.jl  smfish.jl  timelapse.jl
    infer/    summaries.jl  abc_smc.jl  telegraph_likelihood.jl  diagnostics.jl
    generators/ random_grn.jl
    io/       tables.jl  hdf5_anndata.jl  newick.jl
  test/       (unit + statistical tests, see 2.11)
  docs/       Documenter.jl site: tutorials, API, theory notes
  paper/      figure scripts, environment, manuscript sources (see Section 3)
  benchmarks/ timing scripts
```

The 25 current source files collapse into the tree above. The old public API (`Donne`, `ssa`, `tauleap`, `exponential_growth`, NamedTuple reactions) is kept in `compat_v1.jl` as thin deprecated wrappers for one release cycle, so the 2022 workflows still run.

### 2.2 Model specification (`model/`)

* `Species`, `Reaction`, `ReactionModel` types replace the `Donne` mutable struct and its fixed 10-column matrices (which silently cap reactants and products at five and index by column position).
* Propensity kinds, replacing the string matching on reaction names in `rates.jl` (`contains(name, "act")`, `"inhib"`, `"comb_a"`, ...):
  * `MassAction` (0th to 3rd order, with the correct combinatorial factors already in `rates.jl`; bimolecular terms scaled by `1/V`).
  * `HillActivation`, `HillInhibition`, `CombinatorialHill` on **concentrations** `X/V` (this is what the `rate[4] = V` mechanism does today, made explicit).
  * `PromoterSwitch` (telegraph: on/off states as a first-class gene attribute rather than species names containing "on"/"off", which is how `exponential_growth.jl` currently detects them).
  * `VolumeScaled` wrapper (transcription proportional to `V`), and `CopyNumber` (gene copies 1 or 2 depending on cell-cycle position, see 2.4).
* `catalyst_bridge.jl`: `ReactionModel(rn::Catalyst.ReactionSystem)` via a package extension, so Catalyst DSL models can be run in populations.
* Named parameters everywhere; a `ParameterSet` that can be varied per cell (extrinsic noise) and per time (drug effects).

### 2.3 Stochastic kernels (`kernels/`)

Consolidation of the current files:

| Current files | v2.0 kernel | Notes |
|---|---|---|
| `ssa.jl`, `ssa2.jl`, `ssa3.jl`, `ssa_switch.jl`, `ssa_switch2/3/4.jl` | `DirectSSA` | One implementation. Direct method with a species-to-reaction dependency graph (what `update_rates.jl` does) and linear search; optional sorted search. Records at fixed output grid or event-wise. |
| `tauleap.jl`, `tauleapswitch.jl`, `tauleapswitch2.jl`, `non_negative_Poisson_tauleap.jl` | `TauLeap` and `HybridSSATau` | Fixed-step Poisson leaping with **rejection-and-fallback to SSA when any species would go negative** (the "switch" idea), plus the critical-reaction partition of Cao, Gillespie and Petzold. |
| `adaptive_tauleap.jl`, `comp_tau.jl`, `comp_g.jl`, `comp_non_tau.jl`, `compute_L.jl`, `HO_reaction.jl` | `AdaptiveTauLeap` | Step-size selection of Cao et al. 2006 kept, bug-fixed (the current `adaptive_tauleap.jl` computes `dt` but still draws Poisson numbers with `data.tau`, and `non_negative_Poisson_tauleap.jl` calls `rand(pois_rand(...))`, which draws a uniform integer instead of a Poisson variate). |

All kernels: preallocated buffers, explicit `AbstractRNG` argument (reproducible seeds), no `deepcopy` per call, no `Any`-typed fields, `Int` counts throughout, and a common `step!(cell, model, kernel, Δt, rng)` interface used by the population loop.

### 2.4 Cells and populations (`cells/`)

Today `exponential_growth.jl` does: exponential volume growth, division when `V > V_f ≈ 2`, division fraction `0.5 + noise`, binomial partitioning of molecules (`getBinomial`), promoter states copied to daughters, and constant population size by replacing random cells with the new daughters. v2.0 makes each of these an explicit, documented option:

* `Cell`: volume, age, generation, birth time, cell id, parent id, clone/barcode id, copy-number state, molecule counts, promoter states, alive flag.
* Growth: exponential (default), with optional per-cell growth-rate heterogeneity (log-normal), and optional drug-dependent growth inhibition.
* Size control: `Sizer` (current behaviour), `Adder`, `Timer`, each with noise.
* Division: partition fraction `f ~ Normal(0.5, σ)` (current) or Beta; partitioning `Binomial(n, f)` (current) with `BetaBinomial` option for clustered partitioning; promoter states inherited by both daughters.
* Gene replication: copy number 1 -> 2 at a set fraction of the cycle, reset at division (needed to reproduce Sukys and Grima's cell-cycle burst-frequency result and to make the observation layer's cell-cycle confound realistic).
* Death: per-cell hazard from the perturbation layer (2.5); optional constant background death.
* Population control: `ConstantN` (Moran-like replacement, current), `FreeGrowth` (branching process, required for persister and kill-curve simulations), `Logistic` (carrying capacity).
* Multithreading across cells between division events (cells are independent there).

### 2.5 Drug and perturbation layer (`perturb/`)

* Dose schedules `d(t)`: constant, pulses, drug holidays, arbitrary piecewise functions; optional one-compartment pharmacokinetics.
* Effect models per cell (all optional and composable):
  * **State-dependent death hazard** `h(d, p) = h_max · d^m / (EC50^m + d^m) · K^q / (K^q + p^q)`, where `p` is the concentration of a resistance-conferring protein. This is the mechanism in your 2020 MGMT model (survival depends on MGMT level) and in the Iyer et al. interpretation (pre-existing states decide fate).
  * **Growth-rate inhibition** `λ(d)`.
  * **Drug-induced transcription** (dose-dependent increase of `k_on` or `k_tx` for a gene), giving drug-induced rather than pre-existing tolerance (the distinction studied by Corigliano et al. and by the Yanai lab's resistance continuum).
  * **Memory modulation**: a perturbation that increases promoter switching rates (the Harmange et al. "disrupt memory, then treat" concept).
* Gene perturbations for causal benchmark generation: knockdown (factor on `k_tx`), knockout, overexpression, applied to arbitrary subsets of cells, with observational/interventional pairing recorded in metadata.

### 2.6 Lineage layer (`lineage/`)

* Full division tree (parent pointers, division times), clone barcodes assigned at a chosen time (MemorySeq/MeRLin-style), Newick and edge-list export.
* Heritability statistics: mother-daughter and sister-sister correlations, lineage autocorrelation and a fitted memory timescale per gene, population-vs-lineage noise decomposition (following Zhang, Singh et al. 2025).
* Luria-Delbrück fluctuation test utilities: grow clones from single cells, compute clone-level variance of a phenotype against the Poisson expectation, and the MemorySeq-style per-gene clonal variance score. Ground-truth "memory genes" are known from the model (slow promoter switching relative to the cell-cycle time).

### 2.7 Observation layer (`observe/`)

* scRNA-seq: per-cell capture efficiency `β_i ~ Beta`, binomial thinning of true counts, Poisson/multinomial sequencing depth, optional extra zero inflation, optional batch effects, and cell-size coupled library size (this falls out automatically because counts scale with volume, and is exactly the confound that transcriptome-size normalisation papers discuss). Snapshot sampling of an asynchronous population at chosen times. Metadata carried with every cell: true volume, age, cycle position, generation, lineage ids, drug exposure history.
* smFISH mode (no capture loss, per-cell absolute counts) and time-lapse reporter mode (protein reporter with measurement noise, sampled at frame intervals), so the same simulation can be "measured" three ways.
* Optional nascent/mature RNA split for RNA-velocity ground truth (D10).

### 2.8 Inference layer (`infer/`)

* Summary statistics: per-gene moments, Fano factor, zero fraction, histogram or Wasserstein distance, per-cycle-phase statistics, lineage correlations.
* `ABC-SMC` (Toni et al. 2009) with adaptive tolerance schedule, using the population simulator as the forward model; parallel over particles.
* Exact telegraph likelihood (Peccoud-Ycart Beta-Poisson) for reference fits without division, to quantify the bias of ignoring division, partitioning and copy number.
* Diagnostics: posterior predictive checks, simple identifiability profiles.

### 2.9 Generators and IO

* `random_grn`: Erdős-Rényi, scale-free and modular topologies, activation/inhibition assignment, optional master regulators; replaces the `@eval`/`global` code in `network_generator.jl` and `random_network.jl`.
* IO: CSV tables (as now), HDF5 in an AnnData-compatible layout (`X`, `obs`, `var`, `obsm`, `uns`) so scanpy and Python benchmark tools can read outputs directly, and Newick trees.

### 2.10 Engineering

* GitHub Actions CI (Linux/macOS, Julia 1.10 and 1), Documenter.jl docs with tutorials, Aqua.jl quality checks, code coverage.
* Reproducibility: every simulation takes an `rng`; population runs record seeds and parameters in the output metadata.
* Deprecation shims for the v1 API with warnings.
* `Project.toml`: compat bounds for all deps; drop `LsqFit`, `StringDistances`, `DelimitedFiles` from hard dependencies (LsqFit is used only for `growth_estimate`, which becomes an optional utility).

### 2.11 Validation tests baked into `test/`

Statistical tests compare simulated stationary distributions to references with Kolmogorov-Smirnov or chi-square at fixed seeds and large samples:

| Test | Reference |
|---|---|
| Birth-death mRNA is Poisson(k/γ) | textbook |
| Telegraph mRNA is Beta-Poisson | Peccoud and Ycart 1995; Raj et al. 2006 |
| Bursty protein is negative binomial / Gamma limit | Friedman et al. 2006; Shahrezaei and Swain 2008 |
| Second-order mass action propensities | combinatorial `n(n-1)/2` factors already in `rates.jl` |
| DirectSSA vs TauLeap vs AdaptiveTauLeap agree within tolerance | internal |
| DirectSSA vs JumpProcesses.jl Direct on the same model | independent implementation |
| Concentration homeostasis under volume-scaled transcription | Bertaux, Marguerat and Shahrezaei 2018 |
| Partitioning-noise contribution to CV² | Huh and Paulsson 2011 |
| Lineage vs population noise differ as predicted | Zhang, Singh et al. 2025 |
| Apparent burst frequency doubles after gene replication | Sukys and Grima 2025 |
| ABC-SMC recovers telegraph parameters within posterior credible intervals | internal, exact likelihood as reference |
| v1 compatibility: `Donne` + `ssa` reproduce v2 `DirectSSA` statistics | regression |

### 2.12 Performance targets

* 10,000 cells × 50 genes × 1,000 output steps in minutes on a laptop (4 threads), with allocation-free inner loops.
* Benchmark table in the paper: throughput vs cells, genes, and kernel; comparison to AgentBasedModeling.jl on a shared 2-gene growth/division model if it installs in the environment.

### 2.13 Risks and mitigations

* **Scope creep**: the layers are ordered; A and B are required for the paper, C and D can be trimmed to supplementary material.
* **Breaking users**: v1 shim plus a migration guide in the docs.
* **Environment limits** (Section 5): analytic and cross-simulator validation and all synthetic benchmarks run here; the real-data supplements (GEO downloads) run on your machine with scripts I provide.
* **Performance regressions**: benchmark suite in `benchmarks/` run in CI on a small model.

---

## 3. Paper plan

### 3.1 Title options

1. "Biomodelling.jl 2.0: mechanistic simulation of heritable gene-expression states, cell division and drug response for single-cell benchmarking and persister biology"
2. "Simulating drug-tolerant persisters from first principles: a multiscale stochastic framework for growing, dividing, drug-treated cell populations"
3. "A mechanistic ground truth for non-genetic drug tolerance: Biomodelling.jl 2.0"

### 3.2 Core claims (each maps to a results section and a figure)

C1. The engine is exact where exact results exist and agrees with an independent simulator (Fig. 2).
C2. Growth, volume-scaled transcription, division and partitioning generate expression memory, cell-size scaling and cell-cycle effects that reproduce published theory, without any additional assumptions (Fig. 3).
C3. Adding a state-dependent drug-death process yields persister dynamics that reproduce the quantitative signatures reported in recent experimental work: dose-dependent population decay with unchanged single-cell division-time distributions, lineage-correlated fates, unchanged barcode diversity for pre-existing states, and optimal intermittent schedules (Fig. 4).
C4. The simulator provides ground truth that existing single-cell simulators cannot: memory genes, lineage-aware benchmarks, division-confounded GRN inference, and causal perturbation data (Fig. 5).
C5. The same forward model supports parameter inference from snapshot and lineage data, and quantifies the bias of ignoring division (Fig. 6).

### 3.3 Approach and scope

* In scope: the five layers above, synthetic validation, reproduction of the authors' 2020 and 2022 results, benchmark demonstrations, one inference demonstration.
* Out of scope (stated in Discussion): spatial structure, cell-cell signalling, mechanistic pharmacology beyond one-compartment PK, tissue-scale agent-based modelling (PhysiCell territory), and deep-learning inference.

### 3.4 Results sections with simulation designs

**R1. Framework and design.** Schematic of the multiscale model; a 15-line code example; feature comparison table against SERGIO, dyngen, scMultiSim, TedSim, LineageSim, GRouNdGAN, AgentBasedModeling.jl and Biomodelling 1.x. (Fig. 1, Table 1)

**R2. Engine validation.** Simulations: birth-death, telegraph, bursty protein, dimerisation network, a 10-gene random GRN. Outputs: stationary distributions vs analytic curves; DirectSSA vs TauLeap vs AdaptiveTauLeap vs JumpProcesses.jl; runtime scaling. (Fig. 2)

**R3. Emergent memory, size scaling and cell-cycle effects.** Simulations: single-gene telegraph in growing/dividing cells across promoter switching rates spanning the cell-cycle time; with/without partitioning noise; with/without copy-number doubling. Outputs: single-lineage traces; mRNA vs volume scaling and concentration homeostasis; mother-daughter and sister correlations vs switching rate; fitted memory timescale vs cell-cycle time; lineage vs population noise; apparent burst-frequency change across the cycle. Comparison targets: Bertaux 2018, Huh and Paulsson 2011, Zhang and Singh 2025, Sukys and Grima 2025. (Fig. 3)

**R4. Drug response and persisters.** Simulations: free-growth population with a resistance gene R (slow-switching promoter) and state-dependent death; three dose levels; schedules: continuous, holidays of varying length, intermediate doses; memory-modulator pretreatment; MGMT/temozolomide re-implementation. Outputs: kill curves and fractional killing; population decay rate vs dose while intermitotic and death-time distributions barely change (target: the 3-fold decay-rate increase with minor single-cell changes reported by Iyer et al.); sister/cousin fate correlations; barcode diversity before and after drug; schedule optimisation heatmap (target: non-zero release periods and lower-than-maximal doses optimal, Corigliano et al.); fewer resistant colonies after memory-disrupting pretreatment (target: Harmange et al.); stable inheritance of high-MGMT state after treatment (target: Lasri and Sturrock 2020). (Fig. 4)

**R5. Benchmark ground truth.** (a) Memory genes: 200-gene random GRN with a known subset of slow-switching genes; run MemorySeq-style clonal variance and a covariance-eigenspectrum power-law test (the Power-Seek principle) and report precision/recall vs memory timescale. (b) GRN inference: same networks simulated with and without division/volume scaling; PIDC, GENIE3/GRNBoost2, correlation baselines; AUPR on raw counts vs concentration-normalised counts. (c) Regression of the 2022 imputation result with two of the same methods. (d) Perturbation ground truth: observational + knockdown datasets; error of additive and mean baselines vs the true shifted distributions. (Fig. 5)

**R6. Inference.** ABC-SMC on synthetic snapshot data: posterior vs exact-likelihood fit for the telegraph model; bias when the fitted model ignores division and partitioning; recovery of drug parameters (EC50, h_max, switching rates) from population kill curves plus lineage correlations. (Fig. 6)

### 3.5 Validation matrix

| Claim | Validation | Source | Data needed | Runs where |
|---|---|---|---|---|
| Exactness | KS tests vs analytic distributions | Peccoud-Ycart; Friedman; Shahrezaei-Swain | none | here |
| Independent implementation | Same model in JumpProcesses.jl | SciML | none | here (if installable), else your machine |
| Size scaling, homeostasis | Reproduce predicted relations | Bertaux et al. 2018 | none | here |
| Partitioning noise | CV² decomposition | Huh and Paulsson 2011 | none | here |
| Lineage vs population noise | Reproduce analytic predictions | Zhang, Singh et al. 2025 | none | here |
| Cell-cycle burst frequency | Reproduce halving per allele after replication | Sukys and Grima 2025 (published numbers) | none; real mESC data optional | here; optional real data on your machine |
| Persister kinetics | Decay rate vs dose; fate correlations; barcode diversity | Iyer, Granada, Chakrabarti 2025 (published numbers) | none | here |
| Schedules | Optimal release period and dose | Corigliano et al. 2025 | none | here |
| Memory disruption | Fewer resistant colonies after pretreatment | Harmange et al. 2023 | none for qualitative; MemorySeq data optional | here; optional on your machine |
| Authors' prior results | MGMT phenotypic selection; 2022 imputation benchmark | Lasri and Sturrock 2020; Lasri, Shahrezaei, Sturrock 2022 | none | here |
| Realism vs other simulators | Distributional metrics vs SERGIO output | SERGIO (git clone works) | none | here |

### 3.6 Figures

* Fig. 1: framework schematic (cell, reactions, growth, division, partitioning, lineage, drug, observation, inference) and code snippet.
* Fig. 2: (a-c) analytic distribution overlays, (d) kernel agreement, (e) JumpProcesses cross-check, (f) runtime scaling.
* Fig. 3: (a) lineage traces of mRNA and volume, (b) count vs volume and concentration homeostasis, (c) mother-daughter and sister correlations vs switching rate, (d) memory timescale vs cycle time, (e) lineage vs population noise, (f) burst frequency across the cycle.
* Fig. 4: (a) schedules and population trajectories, (b) kill curves, (c) decay rate vs dose with division-time distributions inset, (d) fate correlations by relatedness, (e) barcode diversity, (f) schedule optimisation heatmap, (g) memory-disruption pretreatment, (h) MGMT reproduction.
* Fig. 5: (a) memory-gene recovery precision/recall, (b) GRN AUPR with/without division and normalisation, (c) 2022 imputation regression, (d) perturbation baseline error.
* Fig. 6: (a) ABC posterior vs exact likelihood, (b) bias from ignoring division, (c) drug-parameter recovery, (d) posterior predictive checks.
* Supplementary figures: algorithm flowcharts, tolerance schedules and convergence, parameter sensitivity, additional kernels, extra runtime tables, API walkthrough.
* Tables: T1 feature comparison; T2 validation matrix; T3 all model parameters with sources; T4 runtime benchmarks.

### 3.7 Methods outline

Model equations (propensities, volume scaling, Hill regulation, promoter switching); growth, size control, division and partitioning rules; copy-number model; population control; drug schedule, PK and effect models; lineage recording and statistics; observation models; kernels and step-size selection; ABC-SMC; benchmark protocols (networks, methods, metrics); parameter tables; software and reproducibility.

### 3.8 Supplement outline

S1 derivations used for validation; S2 kernel details and correctness tests; S3 full parameter tables; S4 additional persister simulations (drug-induced vs pre-existing tolerance); S5 benchmark details; S6 inference diagnostics; S7 runtime; S8 migration guide from v1.

### 3.9 Value and novelty statement

Existing simulators either model regulation without growth and division (SERGIO, dyngen, scMultiSim, GRouNdGAN) or model lineage trees phenomenologically without mechanism (TedSim, LineageSim). AgentBasedModeling.jl provides mechanistic growth and division but no drug layer, lineage statistics, observation model, inference or benchmark tooling. Biomodelling.jl 2.0 is the only tool combining all five in one place, targeted at the problem the field is actively working on: non-genetic drug tolerance driven by heritable expression states.

### 3.10 Reproducibility package

`paper/` with one script per figure, a `Project.toml`/`Manifest.toml` pin, a Python `requirements.txt` for benchmark tools, fixed seeds, a `make figures` entry point, and a Zenodo DOI for the release.

### 3.11 bioRxiv submission checklist

Category Systems Biology; CC-BY license; code at the GitHub release and Zenodo DOI; data availability statement (all data synthetic and regenerable, real-data supplements cite their GEO accessions); author contributions; competing interests; funding.

---

## 4. Execution plan and review checkpoints

| Phase | Deliverable | Review checkpoint |
|---|---|---|
| A. Foundation | Restructure, model types, kernels, cells layer, tests, CI, docs skeleton, v1 shim, LICENSE | Checkpoint 1: you review the API and test results |
| B. Biology layers | Drug/perturbation, lineage, observation, IO | Checkpoint 2: you review a persister demo notebook |
| C. Inference and generators | ABC-SMC, telegraph likelihood, random GRNs, benchmark scripts | Checkpoint 3 |
| D. Paper simulations | `paper/` figure scripts, all six figures, draft manuscript (Markdown + BibTeX, exportable to PDF/DOCX), supplement | Checkpoint 4: full draft for your edits |
| E. Finalisation | Revisions, optional real-data supplements run on your machine, release tag, Zenodo, submission package | Submission |

Each phase ends with a commit on the branch and a summary of what changed and what was verified.

---

## 5. Environment notes (verified in this session)

* Julia's official download hosts and package server are blocked by the network policy; conda-forge is reachable, so Julia is being installed through it. Package sources are reachable through git; binary artifacts may not be. Pure-Julia dependencies are therefore preferred throughout.
* Publisher sites, PubMed/PMC, GEO and Zenodo are blocked from this environment; PyPI and git are reachable. Real-data supplements will be scripted for your machine.
* Literature research for this plan was done with the bioRxiv connector and web search; the papers cited above were read at the abstract level plus published summaries, so exact numbers used as targets will be re-checked against the full texts on your side before they go into the manuscript.
