---
title: "Supplementary Information"
author: "Ayoub Lasri"
date: "Berache Limited, Dublin, Ireland"
---

**Mechanistic simulation of heritable expression states, cell division and drug response in single-cell populations.** Supplementary Notes 1–15, Supplementary Figures 1–6 and Supplementary Tables 1–5.

# Supplementary Note 1: Stationary laws used for validation

**Birth-death.** For $\emptyset \xrightarrow{k} X \xrightarrow{\gamma} \emptyset$ the stationary distribution is Poisson with mean $k/\gamma$.

**Telegraph model.** With activation $k_{\mathrm{on}}$, inactivation $k_{\mathrm{off}}$, transcription $k_{\mathrm{tx}}$ from the active state and degradation $\gamma$, and $a = k_{\mathrm{on}}/\gamma$, $b = k_{\mathrm{off}}/\gamma$, $\lambda = k_{\mathrm{tx}}/\gamma$, the stationary distribution is the Beta-Poisson mixture
$$P(n) = \frac{\lambda^n}{n!}\frac{\Gamma(a+n)\,\Gamma(a+b)}{\Gamma(a+b+n)\,\Gamma(a)}\,{}_1F_1(a+n;\,a+b+n;\,-\lambda),$$
whose mean is $\lambda a/(a+b)$ and whose variance is $\mu + \lambda^2 ab/[(a+b)^2(a+b+1)]$. We evaluate ${}_1F_1(\alpha;\beta;-\lambda)$ through Kummer's transformation $e^{-\lambda}{}_1F_1(\beta-\alpha;\beta;\lambda)$, whose series has positive terms; partial sums are rescaled by $10^{-200}$ whenever they exceed $10^{200}$ so that the log of the sum is accumulated without overflow.

**Bursty protein.** Protein produced in geometric bursts of mean size $b$ arriving at rate $a$ and degraded at rate $\gamma$ is negative binomial with shape $a/\gamma$ and success probability $1/(1+b)$. The package emulates bursts through an mRNA of lifetime $1/(200\gamma)$ translated at rate $200\gamma b$.

# Supplementary Note 2: Kernel details and correctness tests

The direct method draws the waiting time from $\mathrm{Exp}(a_0)$ and the channel by linear search over the cumulative propensities; after each event only the propensities of reactions that depend on a changed species (the dependency graph) are recomputed, and the running total is refreshed every 512 events. Fixed-step tau-leaping draws $\mathrm{Poisson}(a_j\tau)$ firings for every channel; the strict variant errors when a species would become negative, the hybrid variant discards the leap and simulates that interval with the direct method. Discarding is conditioned on the realised leap, so the accepted steps of the hybrid variant follow a truncated Poisson law and the kernel is approximate; Fig. 2d measures the resulting discrepancy against the exact stationary laws at two step sizes. Adaptive tau-leaping classifies as critical every channel within $n_c = 10$ firings of exhausting a reactant, selects $\tau$ from the non-critical channels through the bounds $\max(\varepsilon x_i/g_i, 1)/|\mu_i|$ and $\max(\varepsilon x_i/g_i, 1)^2/\sigma_i^2$, draws the time to the next critical event from $\mathrm{Exp}(a_0^{\mathrm{crit}})$, fires at most one critical channel, halves $\tau$ on rejection, and performs 100 exact events whenever $\tau < 10/a_0$.

Supplementary Table 1 lists the correctness tests of the package test suite. They run on every commit under continuous integration, so the tolerances below are enforced rather than reported once.

Supplementary Table 1. Correctness tests of the kernels and of the population layer (`test/test_kernels.jl`, `test/test_population.jl`, `test/test_exact_population.jl`).

| Model | Kernel | Cells | Statistic | Tolerance |
|---|---|---|---|---|
| Birth-death, $k = 12$, $\gamma = 1$ | direct | 8,000 | Kolmogorov-Smirnov distance to $\mathrm{Poisson}(12)$ | < 0.015 |
| same | hybrid ($\tau = 0.02$), adaptive, strict tau-leap ($\tau = 0.005$) | 8,000 each | same | < 0.02 |
| same | all four | 8,000 each | error of the mean; minimum count | < 0.2; $\geq 0$ |
| Telegraph, $k_{\mathrm{on}} = 0.4$, $k_{\mathrm{off}} = 0.6$, $k_{\mathrm{tx}} = 12$ | direct | 12,000 | Kolmogorov-Smirnov distance to the Beta-Poisson law | < 0.015 |
| same | direct | 12,000 | error of the mean transcript number; of the active-promoter fraction | < 0.15; < 0.02 |
| Bursty protein, $a = 1.5$, $b = 6$ | direct | 8,000 | Kolmogorov-Smirnov distance to $\mathrm{NegBin}(a/\gamma, 1/(1+b))$ | < 0.03 |
| Stiff dimerisation ($k = 50$, $k_d = 0.05$, $k_u = 1$, $g = 0.5$) | adaptive ($\varepsilon = 0.02$) against direct | 2,000 each | relative difference of the mean dimer number | < 5% |
| Birth-death | strict tau-leap ($\tau = 0.5$) on a model that would go negative | 400 | an error is raised rather than a negative count returned | raised |
| Birth-death | direct, one thread against several | 300 | identical final states | exact |
| Telegraph in a dividing population | hybrid | 200 | population size held constant under `ConstantN`; doubling time under `FreeGrowth` | exact; within 10% |
| Constitutive production with a memoryless timer | direct | 3,000 lineages | mean and standard deviation against the closed-form single-lineage law | < 0.5 each |
| same | direct | branching population | mean against the closed-form population law; separation from the lineage law | < 1.0; > 1.5 |
| Gamma interdivision times, $\mathrm{cv} = 0.25$ and $1$ | — | 20,000 draws | error of the mean and of the coefficient of variation | < 0.15; < 0.05 |

# Supplementary Note 3: Parameters of all simulations

Supplementary Table 2. Parameters of every simulation in this Article.

| Figure | Model | Parameters |
|---|---|---|
| 2a | birth-death | $k = 12$, $\gamma = 1$, $T = 15$ |
| 2b | telegraph | $(k_{\mathrm{on}}, k_{\mathrm{off}}, k_{\mathrm{tx}}) = (0.4, 0.6, 12)$ and $(0.05, 0.5, 40)$, $\gamma = 1$, $T = 40$ |
| 2c | bursty protein | $a = 1.5$, $b = 6$, $\gamma = 1$, $T = 25$ |
| 2d-e | telegraph | $(0.4, 0.6, 12)$, $n = 20\,000$ |
| 2f | random telegraph GRNs | $n_{\mathrm{act}} = n_{\mathrm{genes}}$, $n_{\mathrm{inh}} = n_{\mathrm{genes}}/2$, $\lambda = \ln 2/20$, $dt = 0.1$, $T = 20$ |
| 3 | telegraph gene with protein | $k_{\mathrm{tx}} = 30$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 4$, $k_{\mathrm{dp}} = 0.2$, $k_{\mathrm{on}} = k_{\mathrm{off}}$ as indicated; $\lambda = \ln 2/20$, sizer $V_{\mathrm{div}} = 2$ (cv 0.05), $\sigma_f = 0.02$, replication at 50% of the cycle where indicated; panel 3d varies the cell-cycle time over 10, 20 and 40; panel 3f uses $(k_{\mathrm{on}}, k_{\mathrm{off}}, k_{\mathrm{tx}}) = (0.5, 1.5, 40)$ with 3,000 cells |
| 3e | stable protein | $k_{\mathrm{on}} = k_{\mathrm{off}} = 5$, $k_{\mathrm{tx}} = 20$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 2$, $k_{\mathrm{dp}} = 0.05$ |
| 4 | resistance gene | as Fig. 3 with $(k_{\mathrm{on}}, k_{\mathrm{off}}) = (0.002, 0.02)$ (memory, $p_{\mathrm{on}} \approx 9\%$) or $(0.05, 0.5)$ (fast, same $p_{\mathrm{on}}$); founders burnt in for 300 time units at constant size; death $h_{\max} = 0.15$, $\mathrm{EC}_{50} = 0.5$, $m = 2$, $K = 150$, $q = 4$; drug from $t = 40$ of the treatment stage; drug-induced variant $(0.0002, 0.02)$ with $k_{\mathrm{on}} \to k_{\mathrm{on}}(1 + 100 d)$; fitness cost variant: growth $\times (1 - 0.5\, c_P^4/(150^4 + c_P^4))$; memory disruption: switching $\times 20$ on $[0, 40)$ or $[0, 120)$; schedule scan: 300 founders, 240 time units, 20-time-unit exposures, five independent seeds per point; cycle-gated variant (Supplementary Fig. 5) adds $\beta = 0.25$, $\varphi_0 = 0.5$, $w = 0.15$ to the death hazard, and is repeated with $h_{\max}$ raised to 0.295 so that the cycle-averaged hazard matches the cycle-blind one |
| 5a | 40 independent telegraph genes | $k_{\mathrm{on}} = k_{\mathrm{off}} \in [10^{-3}, 1]$ (log-spaced), $k_{\mathrm{tx}} = 20$, $k_{\mathrm{dm}} = 0.5$; 40 founders, 6 doublings |
| 5b-c | random GRN | 30 genes, 30 activations, 15 inhibitions, $k_{\mathrm{tx}} \in [5, 30]$, $K \in [2, 10]$, $n = 2$, basal 0.05; 1,500 cells; sequencing capture 0.15 (cv 0.3) |
| 5d | random GRN | 20 genes, 22 activations, 10 inhibitions; knockdown factor 0.05; 800 cells |
| 6a-b | telegraph | truth $(0.3, 0.6, 20)$; priors log-uniform on $[0.01, 10]$, $[0.01, 10]$, $[1, 200]$; 200 (a) / 150 (b) particles |
| 6c | resistance gene | truth $h_{\max} = 0.5$, $K = 150$; priors log-uniform on $[0.05, 5]$ and $[20, 1000]$; 100 particles |
| 7 | resistance gene, hours | cell cycle 24 h (sizer, cv 0.1), $k_{\mathrm{tx}} = 30$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 4$, $k_{\mathrm{dp}} = 0.2$ per h, $K = 150$, $q = 4$, growth-arrest Hill coefficient $m_g = 2$; fitted: $k_{\mathrm{on}}$, $k_{\mathrm{off}} \in [10^{-3}, 0.3]$, $h_{\max} \in [0.005, 0.5]$ per h, $\mathrm{EC}_{50} \in [3, 40]$ µM, death Hill coefficient $m \in [1, 6]$, $\mathrm{IC}_{50} \in [1, 40]$ µM; 300 founders burnt in for 10 cycles, 48 h drug-free, 72 h drug, integration step $dt = 0.5$ h; 160 Latin-hypercube points and 80 Nelder-Mead iterations with seed 1; identifiability (Supplementary Fig. 3) by 11-point log-grid conditional profiles of the six parameters, a $7\times7$ memory × resistant-fraction slice, and the step halved to 0.25 h with three seeds, via the `profile` part of the script; the cycle-dependent variant (Supplementary Fig. 4) adds a seventh parameter $\beta \in [0, 1]$ with $\varphi_0 = 0.5$, $w = 0.15$, fitted under the same budget by the `cycle` part |
| 8a-e | melanoma-like | $k_{\mathrm{tx}} = 3$, $k_{\mathrm{dm}} = 0.1$, $k_{\mathrm{tl}} = 0.4$, $k_{\mathrm{dp}} = 0.02$ per h; net doubling 4 weeks, memory 5 net doublings, pre-resistant fraction 0.005 (1:200, above the traced 1:1,000 to 1:10,000 of @emert2021; Methods), $h_{\max} = 0.0023$ per h, $\mathrm{EC}_{50} = 0.3$, $m = 2$, protection $K = 150$, $q = 4$, growth arrest $\mathrm{IC}_{50} = 0.3$, $m_g = 2$, with the same protection, $k_{\mathrm{off}} \to k_{\mathrm{off}}/(1 + 9d)$, fitness cost 0.5, partial protection: additional unprotected growth inhibition with $\mathrm{IC}_{50} = 1$; 2,000 founders with promoter states drawn from the stationary distribution, 60 weeks, dt 4 h; progression from the nadir at 1.73 × the running minimum after the 8-week lead-in, loss of control at 1.2 × the pre-treatment size; cycle-gated variant (Supplementary Fig. 5) adds $\beta = 0.25$, $\varphi_0 = 0.5$, $w = 0.15$, and is repeated with $h_{\max} = 0.0045$ per h so that the cycle-averaged hazard matches the cycle-blind one; the memory of the resistant state is scanned over 2, 3.5, 5, 8 and 12 net doublings with three seeds per point |
| S6 | exact-solution validation | constitutive $k = 20$, $\gamma = 1$; bursts of $b = 4$ at $k/b$; telegraph $(k_{\mathrm{on}}, k_{\mathrm{off}}, k_{\mathrm{tx}}) = (0.5, 1, 40)$, all with a memoryless interdivision time of mean 1, $dt = 0.002$, exact kernel, 40,000 lineages and five populations capped at 15,000 cells, burnt in for 12 generations; replication model $k = 20$, $\gamma = 1$, volume-scaled synthesis, exponential growth $\ln 2$, deterministic cycle of 1, replication at 50% of the cycle, founders placed uniformly on the update grid, 14 generations; step-size scan over $dt = 0.02, 0.01, 0.005, 0.002$ |
| 8f-h | MGMT model | same expression kinetics; net doubling 40 days, memory 4 net doublings, MGMT-expressing fraction 0.01 or 0.30, $h_{\max} = 0.03$ per h at peak, $\mathrm{EC}_{50} = 0.4$ of the standard bolus peak, $m = 2$, $K = 150$, $q = 4$; elimination half-life 2.1 h; stoichiometric MGMT consumption at rate $k\,d\,V\,c/(c + K_m)$ with $k = 300$, $K_m = 150$, $k$ chosen so that the consumption rate at the standard bolus peak matches the first-order parameterisation it replaces rather than fitted to data, and with no lesion pool as a state variable, so lesions are assumed to form in proportion to dose and to be repaired at once; 1,500 founders at the stationary promoter distribution, six 28-day cycles, dt 1 h |

# Supplementary Note 4: Additional persister simulations

![](figS1_timing.png)

**Supplementary Fig. 1 | Single-cell timing distributions underlying Fig. 4c.** **a**, Division times of cells born under drug and **b**, times from drug start (or from birth, for cells born under drug) to death, per dose.

The division-time distribution is set by the sizer and is the same at every dose, because in this model the drug kills but does not slow growth. Times to death are broad (coefficient of variation {{death_cv_0_5}} at dose 0.5 and {{death_cv_2_0}} at dose 2) and their mean shifts by much less than the population decay rate, because most deaths occur among low-expressing cells whose hazard is already close to $h_{\max}$ at dose 0.5 ($\mathrm{EC}_{50} = 0.5$, $m = 2$); raising the dose mainly shortens the survival of cells with intermediate protection. The dose dependence of population decay is therefore carried by the fraction of cells that are protected, not by the kinetics of death of unprotected cells, which is the interpretation @iyer2025 give of their single-cell tracking data.

Under continuous dosing at the highest dose (release period 0, dose 2 in Fig. 4h), the long-term growth rate is {{cont_pre}} per time unit for pre-existing tolerance, {{cont_cost}} when the resistant state carries a 50% growth cost and {{cont_ind}} for drug-induced tolerance; the schedules with the lowest long-term growth rate are {{best_pre}}, {{best_cost}} and {{best_ind}} respectively. For pre-existing tolerance without a cost, every release period raises the net growth rate; with the fitness cost, release periods of 5 and 10 time units are within 0.003 per time unit of continuous dosing at doses of 1 and 2; for drug-induced tolerance, dose 1 gives a lower net growth rate than dose 2 at every release period, and release periods of 20 time units or more raise the growth rate in every model. Each point is a population of 300 founder cells. {{sched_seed_note}}

# Supplementary Note 5: Benchmark details

Network inference methods: absolute Pearson and Spearman correlations of $\log(1+x)$ (Pearson) or ranks (Spearman); GENIE3 as random-forest importances (200 trees, `sqrt` features) fitted on $\log(1+x)$ with scikit-learn. Imputation: kNN-smoothing (one step, $k = 15$, 10 principal components of the Freeman-Tukey transformed, library-size normalised counts) and MAGIC with default parameters on square-root normalised counts. Metrics: area under the precision-recall and receiver operating characteristic (ROC) curves against the undirected (correlation) or directed (GENIE3) ground-truth adjacency. The cell-cycle-regressed dataset is the residual of $\log(1+c)$ (concentration $c$) after ordinary least-squares regression on gene copy number, cell age and age squared, per gene.

# Supplementary Note 6: Inference diagnostics

**Summary statistics and distances.** For Fig. 6a,b the distance between a simulated and the observed sample of counts is their Wasserstein-1 distance [@bernton2019], computed as the mean absolute difference between the two empirical quantile functions on a grid of 201 probabilities, divided by the observed mean; each simulated sample contains 600 cells (a) or the final snapshot of 500 founders grown for 80 time units (b), against 2,000 observed cells. The distance between two independent samples at the true parameters (the noise floor) is reported in the summary tables of the repository. For Fig. 6c the summaries of a treated population are the log surviving fraction on a grid of 25 time points from the start of the drug to 60 time units later (populations that went extinct contribute a floor of $10^{-3}$) and the concordance of death fates between sisters born in the 20 time units before the drug; the distance is the root-mean-square difference of these summaries. Each generation's tolerance is the median of the previous generation's accepted distances, particles are perturbed with a Gaussian kernel of variance twice the weighted variance of the previous generation, and the first generation samples the prior.

![](figS2_abc_schedules.png)

**Supplementary Fig. 2 | ABC-SMC diagnostics.** **a**, Tolerance and **b**, acceptance rate per generation for the three runs of Fig. 6.

The run on non-dividing cells ({{fig6a_gens}} generations, 200 particles) reached a final tolerance of {{fig6a_eps}} with an acceptance rate of {{fig6a_acc}} in the last generation; the division-aware run ({{fig6b_gens}} generations, 150 particles, each particle simulating a population of 500 founders for 80 time units) reached {{fig6b_eps}} at {{fig6b_acc}}; the drug-parameter run ({{fig6c_gens}} generations, 100 particles, 200 founders each) reached {{fig6c_eps}} at {{fig6c_acc}}. The importance weights of the final generation are uneven, so fewer particles carry the posterior than were run: by Kish's formula the effective sample sizes of the final generations are {{fig6a_ess}} for the non-dividing run, {{fig6b_ess}} for the division-aware run and {{fig6c_ess}} for the drug-parameter run. The posteriors of Fig. 6a,b are therefore drawn with a few tens of effectively independent particles, and the isolated tall bars separated by empty bins in those histograms are the weight of one or two particles rather than structure in the posterior. They are not a boundary effect: the accepted values run over {{fig6a_span}} in the non-dividing run and {{fig6b_span}} in the division-aware run, against prior supports of 0.01 to 10 on both switching rates and 1 to 200 on the transcription rate, so no particle is clipped at a prior edge and the histograms of Fig. 6a,b are binned over the whole prior support. In the drug-parameter run the accepted values are {{fig6c_span}}, and the upper end of $h_{\max}$ does approach the prior bound of 5, because the ridge of Fig. 6c runs to the edge of the prior. Posterior predictive checks were made by simulating the forward model at the posterior median (Fig. 6d for the drug parameters); the posterior median rather than the mean is reported throughout because the posteriors are asymmetric on the log scale.

# Supplementary Note 7: Runtime

Supplementary Table 3 lists the wall-clock times of the population simulations of Fig. 2f (200 steps of $dt = 0.1$ for random telegraph networks of 10, 50 and 100 genes, four threads of an Intel Xeon processor at 2.1 GHz, Julia 1.10). Throughput in cell-steps per second falls with the number of reactions because every step of every cell has to sample the events of all channels; the hybrid kernel is 2 to 5 times faster than the direct method at these tolerances, the gap widening with the number of reactions. The single-cell interface used for the comparison with JumpProcesses.jl (Fig. 2e) allocates a workspace per simulation and is slower in serial execution than JumpProcesses.jl on that one-gene model; the population interface reuses one workspace per thread.

Supplementary Table 3. Runtime of population simulations.

{{runtime_table}}

# Supplementary Note 8: Laboratory and clinical reference data

Supplementary Table 4. Reference values used for calibration and validation (transcribed from the cited articles; files in `paper/data/`).

| Quantity | Value | Source |
|---|---|---|
| U2OS cells tracked at cisplatin addition (7 / 10 / 13 µM) | 232 / 240 / 296 | @iyer2025, Table 2 |
| of which died within 72 h | 30 / 66 / 176 | same |
| of which divided within 72 h | 100 / 75 / 67 | same |
| of which survived without dividing | 102 / 99 / 53 | same |
| HCT116 cells at drug addition; death fraction | 275; 0.64 | @iyer2025, Table 1 |
| Plateau of the HCT116 kill curve | after ~100 h; the U2OS curves show no biphasic decay at any of the three concentrations | @iyer2025, Fig 1b,e |
| Lineage correlations of fate (HCT116 only; the U2OS dataset tracked one randomly chosen daughter per division, so sister and cousin pairs cannot be formed in it) | present for first and second cousins; absent for third cousins (reported graphically; no numeric coefficients given); the memory of at least two to three generations these imply sets the prior used in the calibration | @iyer2025, Fig 5a |
| Single-cell intermitotic and apoptosis times across doses (U2OS) | not significantly different (Kruskal-Wallis $P = 0.22$ and $P = 0.53$) while population decay rates differ about threefold | @iyer2025 |
| Pre-resistant melanoma cells traced back from resistant fates | initial frequency ~1:1,000 to 1:10,000 (context, not a value used: the melanoma case study uses 1:200; Methods) | @emert2021 |
| Pre-resistant melanoma cells | 1:50 to 1:500 per marker; EGFR-high cells give 7.9 ± 0.9 fold more resistant colonies | @shaffer2017 |
| N15-0385 glioblastoma doubling time | 50 h | @lasri2020 |
| Temozolomide elimination half-life in plasma | 1.8 h (single dose); 2.1 h (population model of plasma and cerebrospinal fluid; the value used here) | @rudek2004; @ostermann2004 |
| Temozolomide penetration of the cerebrospinal fluid | exposure 20% of plasma exposure | @ostermann2004 |
| MGMT / alkyltransferase depletion in peripheral blood mononuclear cells | $-63\%$ at 14 days, $-73\%$ at 21 days on protracted schedules; nadir 18.0 ± 2.26% of initial on a compressed 1,000 mg/m² schedule | @tolcher2003; @middleton2000 |
| Tumour MGMT activity in orthotopic GBM43 xenografts | depleted by day 6 on both schedules; still suppressed at day 22 only on the 21-day schedule; back to baseline in both by day 29 | @robinson2010 |
| RTOG 0525 regimens | 150-200 mg/m² days 1-5 vs 75-100 mg/m² days 1-21 of 28-day cycles; median OS 16.6 vs 14.9 months | @gilbert2013 |
| SWOG S1320 regimens | continuous vs 3 weeks off / 5 weeks on after an 8-week lead-in; median PFS 9.0 vs 5.5 months (HR 1.36 intermittent:continuous, $P = 0.063$; the trial pre-specified two-sided $\alpha = 0.2$ and 80% confidence intervals); median OS 29.2 months in both arms, a secondary end point the trial was not powered for | @algazi2020 |

# Supplementary Note 9: Migration from version 1

The v1 API (`Donne`, `ssa`, `tauleap`, `tauleapswitch`, `adaptive_tauleap`, `exponential_growth`, NamedTuple reactions) remains available through a deprecated compatibility layer that converts NamedTuple reactions to `ReactionModel`s (reactions whose names contain `act`, `inhib`, `comb_a` or `comb_i` become Hill kinetics with `rate = [k, n, K]`; species whose names contain `on`/`off` become promoter groups) and maps the algorithms onto the new kernels. Output shapes are preserved, including the leading `:NULL` column. See `docs/src/migration.md`.

# Supplementary Note 10: Relation to previous work

The ingredients of this framework have precedents, and several of its results confirm conclusions reached earlier by other means. The table separates, claim by claim, what was already established from what this work adds. Prior results the framework reproduces are evidence that it is behaving correctly, not claims of novelty.

| Result in this Article | Already established | What this work adds |
|---|---|---|
| Stochastic reaction kinetics inside growing, dividing cells (Fig. 3) | Simulators of stochastic expression with growth and division exist [@bertaux2018; @piho2025abm; @thomas2021], including version 1 of this framework [@lasri2022] | Gene replication, promoter inheritance, a drug layer, lineage statistics, observation models and inference in one forward model, so that a resistance phenotype is a heritable expression state with measurable kinetics rather than an assumed compartment |
| Synthetic scRNA-seq from a known network in a dividing population, for imputation and network-inference benchmarking (Fig. 5) | Version 1 of this framework already did this [@lasri2022] | Gene replication and cell-cycle copy number, sizer and adder division rules in place of a timer, inherited promoter states, the lineage table, the drug layer, likelihood-free inference and the schedule optimiser; nothing in Fig. 5 that version 1 could also produce is presented as new |
| Exact stationary distributions for growing, dividing cells differ between a single lineage and a population snapshot (Supplementary Note 14) | Solved exactly for bursty expression with replication, partitioning and general interdivision times [@beentjes2020], for the extended telegraph model with volume-dependent synthesis and size control [@jia2023], and framed as an ergodic principle [@thomas2017; @thomas2021] | The population and lineage layers are validated against those solutions in both modes, rather than only against single-cell stationary laws |
| Inference of kinetic parameters from lineage-resolved data | Finite-state-projection inference from mother machines and lineage trees [@piho2024feedback] | Inference of drug parameters from population data, with the identifiability limit made explicit (Fig. 6c) |
| Simulation of division trees for lineage-tracing methods | TedSim couples expression to division history [@pan2022tedsim]; Cassiopeia simulates topologies, heritable fitness and CRISPR barcodes [@jones2020cassiopeia] | Molecular content on the same tree: partitioning at division, promoter states, drug-induced death, and per-cell fates with the molecules that caused them |
| Heritable expression states decide which cells survive a drug (Fig. 4) | Measured directly in barcoded and time-lapse experiments [@shaffer2017; @harmange2023; @iyer2025; @oren2021] | A generative model in which memory is one promoter timescale, reproducing the reported signatures quantitatively and predicting which of them discriminate pre-existing from induced tolerance |
| Cell-cycle-aware inference of transcriptional kinetics is necessary (Fig. 3f, Fig. 6b) | Established from data and theory [@sukys2025; @zhang2025; @okochi2026scdivide] | A simulator that generates the data such methods assume, and a quantification of the bias incurred by a division-blind fit |
| Benchmarks of network inference from single-cell data | Extensive benchmarks exist [@pratapa2020; @dibaeinia2020; @lasri2022] | The growth and cell-cycle confound with known ground truth: how much accuracy is lost in a dividing population, and how much cell-cycle regression recovers (Fig. 5b) |
| The benefit of treatment holidays depends on a cost of resistance and on turnover (Fig. 8a-e) | Compartment models of adaptive therapy [@zhang2017adaptive; @strobl2021] | The same dependence from single cells whose resistance is a heritable expression state; cost and degree of protection become cell properties that lineage experiments measure; partial protection reverses the ranking even without a cost |
| Perturbing the retention of a resistant state changes the outcome (Fig. 4f) | Phenotypic-switching theory [@gunnarsson2020] and memory-disrupting compounds [@harmange2023] | The timing requirement (disruption must continue during exposure) and a lineage-resolved readout, the number of surviving clones |
| An intermediate dose can be optimal against drug-induced tolerance (Fig. 4h) | Reported for induced persisters [@corigliano2025] | Reproduced; the accompanying benefit of release periods does not appear here, which the text attributes to the decay of the induced state in this model rather than to a disagreement about data |
| MGMT expression selects glioblastoma cells under temozolomide (Fig. 8f-h) | Phenotypic selection with stable inheritance [@lasri2020] | Pharmacokinetics, the RTOG 0525 regimens, consumption of MGMT by the drug, and the fractionation that minimises the final population at fixed cumulative dose |
| Dose-dense temozolomide does not improve survival, in methylated and unmethylated tumours alike | Clinical result [@gilbert2013] | A mechanistic account for unmethylated tumours, where the model reproduces the null with MGMT stable and breaks it only when the drug consumes MGMT, which makes the per-lesion consumption rate the measurement that separates the regimens there; in methylated tumours the model predicts a dose-dense advantage from cumulative dose alone, with MGMT stable, a divergence from the trial that consumption deepens but does not cause and that bounds the model's dose response instead |

# Supplementary Note 11: Identifiability and numerical robustness of the calibration

Six parameters were fitted to six measured fractions, so the fit cannot determine all six. Three checks bound what the calibration does and does not establish.

![](figS3_identifiability.png)

**Supplementary Fig. 3 | What the fate fractions determine.** **a**, Root-mean-square error of the simulated fate fractions at the training concentrations when one parameter is moved away from its fitted value over its full range and the others are held at theirs, plotted against the ratio to the fitted value; the dotted line marks one binomial standard error of the measurement above the best fit. **b**, The same error over the memory of the resistant state and the fraction of cells in it, with the death parameters held at their fitted values; red outlines mark the combinations that lie within one standard error of the best. **c**, Simulated fates at 13 µM with the integration step halved (mean ± s.d. of three seeds).

The death parameters are the better determined: the maximal death rate stays within one standard error of the best fit over about a twofold range, and the EC$_{50}$, Hill coefficient and growth-arrest IC$_{50}$ over narrower ranges still (Supplementary Fig. 3a). The two promoter switching rates are not determined: each can move by roughly sixfold with no penalty the data can detect. The two-dimensional slice shows why (Supplementary Fig. 3b). A short memory with many resistant cells and a long memory with few produce nearly the same fate fractions, so the error surface is flat along that direction, and a majority of the combinations tested (memory one to nine generations, resistant fraction 1% to 40%) fit within the noise of the measurement; those that do span one to nine generations of memory and 1% to 25% of cells resistant. The memory quoted in the main text is therefore set by the prior taken from the lineage correlations, not by the fate fractions, and the resistant fraction that accompanies it should not be read as a measurement.

Halving the integration step from 0.5 h to 0.25 h changes each simulated fate fraction by less than the seed-to-seed spread (Supplementary Fig. 3c), so the discretisation of the death hazard is not a source of error at this step size.

A profile in which the remaining parameters are re-optimised at each fixed memory would be preferable to a slice, but at seven seconds per objective evaluation the inner search could not be given enough budget to be reliable: it returned fits three times worse than the known solution at short memory, which reflects the optimiser rather than the data. The slice avoids that failure at the cost of conditioning on the fitted death parameters.

# Supplementary Note 12: Cell-cycle-dependent killing as an alternative to heritable expression

The lethality of cisplatin is enhanced during replication, so where a cell sits in its division cycle when the drug arrives is a competing explanation for two of the observations the calibrated model reproduces: the near-invariance of the mean time to death across concentrations, and the large fraction of cells that survive three days of drug without dividing. A model in which killing is confined to a window of the cycle would spare exactly those cells that never enter that window, without any heritable expression state.

![](figS4_cycle.png)

**Supplementary Fig. 4 | Cell-cycle-dependent killing tested against the fate data.** **a**, Root-mean-square error of the fate fractions at the training concentrations as the cycle-independent fraction $\beta$ of the death hazard is moved from 0 (killing confined to the replication window) to 1 (cycle-independent), with the remaining parameters held at the cycle-dependent fit; dotted line, one binomial standard error of the measurement above the best fit; dashed line, the fitted value. **b**, Simulated against observed fate fractions at all three concentrations for the cycle-independent and the cycle-dependent fit (squares, the held-out 10 µM condition; line, equality).

{{cycle_note}}

The two mechanisms are not exclusive, and the experiment that separates them is not a fate count: it is a measurement of cell-cycle position at the moment the drug arrives, which time-lapse imaging with a cycle reporter provides directly. Until that is done, the heritable-expression reading of the calibrated parameters should be understood as one of two accounts that these fate fractions cannot tell apart.

# Supplementary Note 13: Cell-cycle-dependent killing in the persister and melanoma case studies

Supplementary Note 12 shows that the cisplatin fate data cannot decide whether cells are spared because they inherited a protective state or because they were outside the replication window. The persister simulations of Fig. 4 and the melanoma simulations of Fig. 8a–e use a hazard with no cell-cycle dependence at all, so the same question applies to them: do their conclusions rest on that choice? We repeated both with three quarters of the death hazard confined to a window around mid-cycle ($\beta = 0.25$, $\varphi_0 = 0.5$, $w = 0.15$), leaving everything else unchanged.

The window position is a modelling choice rather than a claim about either agent. The persister study uses a generic drug and commits to no mechanism; the melanoma study is built on a BRAF/MEK inhibitor, which does not create replication-coupled lesions as a platinum drug does, although it does act on cells that are traversing the cycle. What these runs test is therefore robustness to a strong cell-cycle gate of any kind, not the phase specificity of a particular drug.

![](figS5_cycle_robustness.png)

**Supplementary Fig. 5 | The case studies under cell-cycle-gated killing.** Three arms throughout: the cycle-blind hazard; three quarters of it gated to a mid-cycle window; and the same gate with the maximal hazard raised so that the cycle-averaged hazard matches the cycle-blind one. **a**, Net population decay rate over the first 30 time units of treatment against dose (solid, left axis) and the mean time to death of killed cells (dashed, right axis). At matched mean hazard the dose response returns to the cycle-blind curve, which is what shows that the shallower response under the gate was the lower average hazard rather than the gating. **b**, Concordance of lineage fate between sisters and between cousins, above the value expected if fates were independent, for the memory gene and the fast-switching control; at matched hazard the gate leaves no excess concordance at all in the control, which has no usable memory. **c**, Weeks after randomisation to loss of control for the three melanoma mechanisms and the three clinical schedules (mean ± s.d. of four seeds, censored at the end of follow-up), cycle-blind (solid), cycle-gated (mid) and gated at matched hazard (pale).

{{cycle_case_note}}

{{cycle_matched_note}}

# Supplementary Note 14: Exact stationary laws for growing, dividing populations

The kernels of Supplementary Note 2 are validated against stationary laws for a single cell at fixed volume. That leaves the layer above them untested against theory: growth, gene replication, partitioning at division, and the difference between watching one lineage and taking a snapshot of a population. Exact solutions for precisely this class of model exist. Beentjes, Perez-Carrasco and Grima solve the stationary distribution with bursty production, DNA replication, binomial partitioning at mitosis and Erlang or general interdivision times, and give different closed forms for the single-lineage and the population-snapshot setting [@beentjes2020]; Jia and Grima solve the extended telegraph model with volume-dependent synthesis, dosage compensation, partitioning, interdivision-time variability and cell-size control, again separately for lineages and populations [@jia2023]. The difference between the two settings is the ergodic principle: a snapshot of a growing population over-weights cells that have recently divided, and therefore recently lost half their molecules [@thomas2017; @thomas2021].

Four analytically solvable members of that class were simulated with the population layer and compared with the exact solutions in both modes (`paper/scripts/fig9_validation.jl` for the simulations, `paper/scripts/exact_solutions.py` for the solutions). The single-lineage mode uses the `MotherMachine` population control, which keeps one daughter at random at every division, so the recorded cells are statistically independent lineages; the population mode uses free branching growth with uniform subsampling, which preserves the snapshot distribution.

**Memoryless interdivision times.** For the first three models the interdivision time is exponential, so the number of molecules is a Markov process on its own and both stationary laws are available in closed form. Writing the generating function as $G(u) = \sum_j c_j u^j$ with $u = z - 1$, the factorial moments obey

$$c_j = \frac{\sum_i \beta_i\, c_{j-i}}{\gamma j + \lambda\,(1 - 2^{-j})} \quad\text{(lineage)}, \qquad c_j = \frac{\sum_i \beta_i\, c_{j-i}}{\gamma j + \lambda\,(2 - 2^{1-j})} \quad\text{(population)},$$

where $\gamma$ is the decay rate, $\lambda$ the division rate and $\beta_i$ the coefficient of $u^i$ contributed by production (for bursts of fixed size $b$ at rate $k_b$, $\beta_i = k_b\binom{b}{i}$). For the telegraph model the same expansion gives a two-by-two linear system per order, because the promoter state is inherited rather than partitioned. The mean is $k/(\gamma + \lambda/2)$ along a lineage and $k/(\gamma + \lambda)$ in a population, so the two differ by a third at $\gamma = \lambda$ and the comparison has ample power to tell them apart. Because the simulator advances cells by a fixed step $dt$ and divides them at the end of a step, the interdivision time it realises is geometric on that grid; the full pmf of the scheme is therefore the leading eigenvector of $[(1-p)I + \kappa p B]\,\mathrm{e}^{dt A}$ with $p = 1 - \mathrm{e}^{-dt/T}$, $A$ the reaction generator, $B$ binomial thinning and $\kappa = 1$ for a lineage or $2$ for a population, which is what the simulations are compared against; it agrees with the closed-form moments above as $dt \to 0$.

**Deterministic cycle with gene replication.** The fourth model has a deterministic cycle, exponential volume growth, volume-scaled synthesis and gene replication at mid-cycle. Production is linear in volume and copy number and decay is first-order, so the count distribution at every cycle phase is exactly Poisson, with a mean that satisfies a linear recursion over the update steps and is halved at division. The two modes then differ only in how cycle phases are weighted, uniformly along a lineage and towards younger cells in a population; the population weights were obtained by iterating the age map of the scheme from the founder ages actually used, so nothing is assumed about convergence to a stable age distribution. The measured age distribution is compared with that prediction alongside the counts.

![](figS6_exact.png)

**Supplementary Fig. 6 | The population and lineage layers against exact solutions.** **a**, Simulated and exact stationary distributions for constitutive production with a memoryless interdivision time, in single-lineage and population-snapshot mode. **b**, Kolmogorov-Smirnov distance to the exact law of the matching mode, and to the exact law of the other mode, for every model. The dotted line is the 99% critical value for an independent sample of that size, which is exact for the lineage bars and an under-estimate for the population ones, whose cells share ancestors. **c**, Kolmogorov-Smirnov distance against the update step: the theoretical curve is the distance between the stationary law of the scheme and that of the continuous-time model, which is what the discretisation costs, and the measured one sits at the Monte Carlo floor throughout. **d**, Measured and predicted distribution of cell-cycle phase, for the model with gene replication. Along a lineage the phase is uniform. In the population it is the two-level step that a deterministic cycle produces: the snapshot is taken half a cycle after a whole number of generations, and the cells younger than that are the ones whose ancestors have divided once more, so they are twice as numerous. The position of the step is set by when the snapshot is taken and has nothing to do with the replication point, which is also at mid-cycle. The prediction is obtained by iterating the age map of the scheme from the founder ages used.

{{exact_note}}

{{exact_table}}

# Supplementary Note 15: Sensitivity of the melanoma schedules to the memory of the resistant state

{{memory_scan_note}} {{memory_scan_reconcile}}
