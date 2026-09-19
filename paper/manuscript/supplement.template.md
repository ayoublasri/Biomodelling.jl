---
title: "Supplementary Information"
author: "Ayoub Lasri"
date: ""
---

**Mechanistic simulation of heritable expression states, cell division and drug response in single-cell populations.** Supplementary Notes 1–11, Supplementary Figures 1–3 and Supplementary Tables 1–4.

# Supplementary Note 1: Stationary laws used for validation

**Birth-death.** For $\emptyset \xrightarrow{k} X \xrightarrow{\gamma} \emptyset$ the stationary distribution is Poisson with mean $k/\gamma$.

**Telegraph model.** With activation $k_{\mathrm{on}}$, inactivation $k_{\mathrm{off}}$, transcription $k_{\mathrm{tx}}$ from the active state and degradation $\gamma$, and $a = k_{\mathrm{on}}/\gamma$, $b = k_{\mathrm{off}}/\gamma$, $\lambda = k_{\mathrm{tx}}/\gamma$, the stationary distribution is the Beta-Poisson mixture
$$P(n) = \frac{\lambda^n}{n!}\frac{\Gamma(a+n)\,\Gamma(a+b)}{\Gamma(a+b+n)\,\Gamma(a)}\,{}_1F_1(a+n;\,a+b+n;\,-\lambda),$$
whose mean is $\lambda a/(a+b)$ and whose variance is $\mu + \lambda^2 ab/[(a+b)^2(a+b+1)]$. We evaluate ${}_1F_1(\alpha;\beta;-\lambda)$ through Kummer's transformation $e^{-\lambda}{}_1F_1(\beta-\alpha;\beta;\lambda)$, whose series has positive terms; partial sums are rescaled by $10^{-200}$ whenever they exceed $10^{200}$ so that the log of the sum is accumulated without overflow.

**Bursty protein.** Protein produced in geometric bursts of mean size $b$ arriving at rate $a$ and degraded at rate $\gamma$ is negative binomial with shape $a/\gamma$ and success probability $1/(1+b)$. The package emulates bursts through an mRNA of lifetime $1/(200\gamma)$ translated at rate $200\gamma b$.

# Supplementary Note 2: Kernel details and correctness tests

The direct method draws the waiting time from $\mathrm{Exp}(a_0)$ and the channel by linear search over the cumulative propensities; after each event only the propensities of reactions that depend on a changed species (the dependency graph) are recomputed, and the running total is refreshed every 512 events. Fixed-step tau-leaping draws $\mathrm{Poisson}(a_j\tau)$ firings for every channel; the strict variant errors when a species would become negative, the hybrid variant discards the leap and simulates that interval with the direct method. Discarding is conditioned on the realised leap, so the accepted steps of the hybrid variant follow a truncated Poisson law and the kernel is approximate; Fig. 2d measures the resulting discrepancy against the exact stationary laws at two step sizes. Adaptive tau-leaping classifies as critical every channel within $n_c = 10$ firings of exhausting a reactant, selects $\tau$ from the non-critical channels through the bounds $\max(\varepsilon x_i/g_i, 1)/|\mu_i|$ and $\max(\varepsilon x_i/g_i, 1)^2/\sigma_i^2$, draws the time to the next critical event from $\mathrm{Exp}(a_0^{\mathrm{crit}})$, fires at most one critical channel, halves $\tau$ on rejection, and performs 100 exact events whenever $\tau < 10/a_0$.

Supplementary Table 1 lists the tests of the package test suite that validate the kernels: Kolmogorov-Smirnov distances below 0.015 (direct) or 0.02 (approximate kernels) to the Poisson, Beta-Poisson and negative-binomial laws with 8 000 to 12 000 cells; agreement of the adaptive kernel with the direct method on a stiff dimerisation network within 5% in the mean; identical results with one and with several threads.

# Supplementary Note 3: Parameters of all simulations

Supplementary Table 2. Parameters of every simulation in this Article.

| Figure | Model | Parameters |
|---|---|---|
| 2a | birth-death | $k = 12$, $\gamma = 1$, $T = 15$ |
| 2b | telegraph | $(k_{\mathrm{on}}, k_{\mathrm{off}}, k_{\mathrm{tx}}) = (0.4, 0.6, 12)$ and $(0.05, 0.5, 40)$, $\gamma = 1$, $T = 40$ |
| 2c | bursty protein | $a = 1.5$, $b = 6$, $\gamma = 1$, $T = 25$ |
| 2d-e | telegraph | $(0.4, 0.6, 12)$, $n = 20\,000$ |
| 2f | random telegraph GRNs | $n_{\mathrm{act}} = n_{\mathrm{genes}}$, $n_{\mathrm{inh}} = n_{\mathrm{genes}}/2$, $\lambda = \ln 2/20$, $dt = 0.1$, $T = 20$ |
| 3 | telegraph gene with protein | $k_{\mathrm{tx}} = 30$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 4$, $k_{\mathrm{dp}} = 0.2$, $k_{\mathrm{on}} = k_{\mathrm{off}}$ as indicated; $\lambda = \ln 2/20$, sizer $V_{\mathrm{div}} = 2$ (cv 0.05), $\sigma_f = 0.02$, replication at 50% of the cycle where indicated |
| 3e | stable protein | $k_{\mathrm{on}} = k_{\mathrm{off}} = 5$, $k_{\mathrm{tx}} = 20$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 2$, $k_{\mathrm{dp}} = 0.05$ |
| 4 | resistance gene | as Fig. 3 with $(k_{\mathrm{on}}, k_{\mathrm{off}}) = (0.002, 0.02)$ (memory, $p_{\mathrm{on}} \approx 9\%$) or $(0.05, 0.5)$ (fast, same $p_{\mathrm{on}}$); founders burnt in for 300 time units at constant size; death $h_{\max} = 0.15$, $\mathrm{EC}_{50} = 0.5$, $m = 2$, $K = 150$, $q = 4$; drug from $t = 40$ of the treatment stage; drug-induced variant $(0.0002, 0.02)$ with $k_{\mathrm{on}} \to k_{\mathrm{on}}(1 + 100 d)$; fitness cost variant: growth $\times (1 - 0.5\, c_P^4/(150^4 + c_P^4))$; memory disruption: switching $\times 20$ on $[0, 40)$ or $[0, 120)$; schedule scan: 300 founders, 240 time units, 20-time-unit exposures |
| 5a | 40 independent telegraph genes | $k_{\mathrm{on}} = k_{\mathrm{off}} \in [10^{-3}, 1]$ (log-spaced), $k_{\mathrm{tx}} = 20$, $k_{\mathrm{dm}} = 0.5$; 40 founders, 6 doublings |
| 5b-c | random GRN | 30 genes, 30 activations, 15 inhibitions, $k_{\mathrm{tx}} \in [5, 30]$, $K \in [2, 10]$, $n = 2$, basal 0.05; 1 500 cells; sequencing capture 0.15 (cv 0.3) |
| 5d | random GRN | 20 genes, 22 activations, 10 inhibitions; knockdown factor 0.05; 800 cells |
| 6a-b | telegraph | truth $(0.3, 0.6, 20)$; priors log-uniform on $[0.01, 10]$, $[0.01, 10]$, $[1, 200]$; 200 (a) / 120 (b) particles |
| 6c | resistance gene | truth $h_{\max} = 0.5$, $K = 150$; priors log-uniform on $[0.05, 5]$ and $[20, 1000]$; 100 particles |
| 7 | resistance gene, hours | cell cycle 24 h (sizer, cv 0.1), $k_{\mathrm{tx}} = 30$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 4$, $k_{\mathrm{dp}} = 0.2$ per h, $K = 150$, $q = 4$, growth-arrest Hill coefficient 2; fitted: $k_{\mathrm{on}}$, $k_{\mathrm{off}} \in [10^{-3}, 0.3]$, $h_{\max} \in [0.005, 0.5]$ per h, $\mathrm{EC}_{50} \in [3, 40]$ µM, $m \in [1, 6]$, $\mathrm{IC}_{50} \in [1, 40]$ µM; 300 founders burnt in for 10 cycles, 48 h drug-free, 72 h drug; 160 Latin-hypercube points and 80 Nelder-Mead iterations with seed 1 |
| 8a-d | melanoma-like | $k_{\mathrm{tx}} = 3$, $k_{\mathrm{dm}} = 0.1$, $k_{\mathrm{tl}} = 0.4$, $k_{\mathrm{dp}} = 0.02$ per h; net doubling 4 weeks, memory 5 generations, pre-resistant fraction 0.005, $h_{\max} = 0.0023$ per h, $\mathrm{EC}_{50} = 0.3$, $m = 2$, protection $K = 150$, $q = 4$, growth arrest $\mathrm{IC}_{50} = 0.3$ with the same protection, $k_{\mathrm{off}} \to k_{\mathrm{off}}/(1 + 9d)$, fitness cost 0.5, partial protection: additional unprotected growth inhibition with $\mathrm{IC}_{50} = 1$; 2 000 founders with promoter states drawn from the stationary distribution, 60 weeks, dt 4 h; progression at 1.73 × nadir after the 8-week lead-in |
| 8e-g | MGMT model | same expression kinetics; net doubling 40 days, memory 4 generations, MGMT-expressing fraction 0.01 or 0.30, $h_{\max} = 0.03$ per h at peak, $\mathrm{EC}_{50} = 0.4$ of the standard bolus peak, $m = 2$, $K = 150$, $q = 4$; elimination half-life 2.1 h; MGMT consumption $k_{\mathrm{dp}} \to k_{\mathrm{dp}}(1 + 20d)$; 1 500 founders at the stationary promoter distribution, six 28-day cycles, dt 1 h |

# Supplementary Note 4: Additional persister simulations

![](figS1_timing.png)

**Supplementary Fig. 1 | Single-cell timing distributions underlying Fig. 4c.** **a**, Division times of cells born under drug and **b**, times from drug start (or from birth, for cells born under drug) to death, per dose.

The division-time distribution is set by the sizer and is the same at every dose, because in this model the drug kills but does not slow growth. Times to death are broad (coefficient of variation {{death_cv_0_5}} at dose 0.5 and {{death_cv_2_0}} at dose 2) and their mean shifts by much less than the population decay rate, because most deaths occur among low-expressing cells whose hazard is already close to $h_{\max}$ at dose 0.5 ($\mathrm{EC}_{50} = 0.5$, $m = 2$); raising the dose mainly shortens the survival of cells with intermediate protection. The dose dependence of population decay is therefore carried by the fraction of cells that are protected, not by the kinetics of death of unprotected cells, which is the interpretation @iyer2025 give of their single-cell tracking data.

Under continuous dosing at the highest dose (release period 0, dose 2 in Fig. 4f), the long-term growth rate is {{cont_pre}} per time unit for pre-existing tolerance, {{cont_cost}} when the resistant state carries a 50% growth cost and {{cont_ind}} for drug-induced tolerance; the schedules with the lowest long-term growth rate are {{best_pre}}, {{best_cost}} and {{best_ind}} respectively. For pre-existing tolerance without a cost, every release period raises the net growth rate; with the fitness cost, release periods of 5 and 10 time units are within 0.003 per time unit of continuous dosing at doses of 1 and 2; for drug-induced tolerance, dose 1 gives a lower net growth rate than dose 2 at every release period, and release periods of 20 time units or more raise the growth rate in every model. All values are single realisations with 300 founder cells.

# Supplementary Note 5: Benchmark details

Network inference methods: absolute Pearson and Spearman correlations of $\log(1+x)$ (Pearson) or ranks (Spearman); GENIE3 as random-forest importances (200 trees, `sqrt` features) fitted on $\log(1+x)$ with scikit-learn. Imputation: kNN-smoothing (one step, $k = 15$, 10 principal components of the Freeman-Tukey transformed, library-size normalised counts) and MAGIC with default parameters on square-root normalised counts. Metrics: area under the precision-recall and receiver-operating curves against the undirected (correlation) or directed (GENIE3) ground-truth adjacency. The cell-cycle-regressed dataset is the residual of $\log(1+c)$ (concentration $c$) after ordinary least-squares regression on gene copy number, cell age and age squared, per gene.

# Supplementary Note 6: Inference diagnostics

**Summary statistics and distances.** For Fig. 6a,b the distance between a simulated and the observed sample of counts is their Wasserstein-1 distance [@bernton2019], computed as the mean absolute difference between the two empirical quantile functions on a grid of 201 probabilities, divided by the observed mean; each simulated sample contains 600 cells (a) or the final snapshot of 500 founders grown for 80 time units (b), against 2 000 observed cells. The distance between two independent samples at the true parameters (the noise floor) is reported in the summary tables of the repository. For Fig. 6c the summaries of a treated population are the log surviving fraction on a grid of 25 time points from the start of the drug to 60 time units later (populations that went extinct contribute a floor of $10^{-3}$) and the concordance of death fates between sisters born in the 20 time units before the drug; the distance is the root-mean-square difference of these summaries. Each generation's tolerance is the median of the previous generation's accepted distances, particles are perturbed with a Gaussian kernel of variance twice the weighted variance of the previous generation, and the first generation samples the prior.

![](figS2_abc_schedules.png)

**Supplementary Fig. 2 | ABC-SMC diagnostics.** **a**, Tolerance and **b**, acceptance rate per generation for the three runs of Fig. 6.

The run on non-dividing cells ({{fig6a_gens}} generations, 200 particles) reached a final tolerance of {{fig6a_eps}} with an acceptance rate of {{fig6a_acc}} in the last generation; the division-aware run ({{fig6b_gens}} generations, 120 particles, each particle simulating a population of 250 founders for 80 time units) reached {{fig6b_eps}} at {{fig6b_acc}}; the drug-parameter run ({{fig6c_gens}} generations, 100 particles, 200 founders each) reached {{fig6c_eps}} at {{fig6c_acc}}. Posterior predictive checks were made by simulating the forward model at the posterior median (Fig. 6d for the drug parameters); the posterior median rather than the mean is reported throughout because the posteriors are asymmetric on the log scale.

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
| Plateau of the HCT116 kill curve | after ~100 h | @iyer2025, Fig 1b |
| Lineage correlations of fate | present for sisters, first and second cousins; absent for third cousins | @iyer2025, Fig 5a |
| Pre-resistant melanoma cells | 1:50 to 1:500 per marker; EGFR-high cells give 7.9 ± 0.9 fold more resistant colonies | @shaffer2017 |
| N15-0385 glioblastoma doubling time | 50 h | @lasri2020 |
| Temozolomide elimination half-life | 2.1 h | @ostermann2004 |
| RTOG 0525 regimens | 150-200 mg/m² days 1-5 vs 75-100 mg/m² days 1-21 of 28-day cycles; median OS 16.6 vs 14.9 months | @gilbert2013 |
| SWOG S1320 regimens | continuous vs 3 weeks off / 5 weeks on after an 8-week lead-in; median PFS 9.0 vs 5.5 months | @algazi2020 |

# Supplementary Note 9: Migration from version 1

The v1 API (`Donne`, `ssa`, `tauleap`, `tauleapswitch`, `adaptive_tauleap`, `exponential_growth`, NamedTuple reactions) remains available through a deprecated compatibility layer that converts NamedTuple reactions to `ReactionModel`s (reactions whose names contain `act`, `inhib`, `comb_a` or `comb_i` become Hill kinetics with `rate = [k, n, K]`; species whose names contain `on`/`off` become promoter groups) and maps the algorithms onto the new kernels. Output shapes are preserved, including the leading `:NULL` column. See `docs/src/migration.md`.

# Supplementary Note 10: Relation to previous work

The ingredients of this framework have precedents, and several of its results confirm conclusions reached earlier by other means. The table separates, claim by claim, what was already established from what this work adds. Prior results the framework reproduces are evidence that it is behaving correctly, not claims of novelty.

| Result in this Article | Already established | What this work adds |
|---|---|---|
| Stochastic reaction kinetics inside growing, dividing cells (Fig. 3) | Simulators of stochastic expression with growth and division exist [@bertaux2018; @piho2025abm], including version 1 of this framework [@lasri2022] | Gene replication, promoter inheritance, a drug layer, lineage statistics, observation models and inference in one forward model, so that a resistance phenotype is a heritable expression state with measurable kinetics rather than an assumed compartment |
| Simulation of division trees for lineage-tracing methods | TedSim couples expression to division history [@pan2022tedsim]; Cassiopeia simulates topologies, heritable fitness and CRISPR barcodes [@jones2020cassiopeia] | Molecular content on the same tree: partitioning at division, promoter states, drug-induced death, and per-cell fates with the molecules that caused them |
| Heritable expression states decide which cells survive a drug (Fig. 4) | Measured directly in barcoded and time-lapse experiments [@shaffer2017; @harmange2023; @iyer2025; @oren2021] | A generative model in which memory is one promoter timescale, reproducing the reported signatures quantitatively and predicting which of them discriminate pre-existing from induced tolerance |
| Cell-cycle-aware inference of transcriptional kinetics is necessary (Fig. 3f, Fig. 6b) | Established from data and theory [@sukys2025; @zhang2025; @okochi2026scdivide] | A simulator that generates the data such methods assume, and a quantification of the bias incurred by a division-blind fit |
| Benchmarks of network inference from single-cell data | Extensive benchmarks exist [@pratapa2020; @dibaeinia2020; @lasri2022] | The growth and cell-cycle confound with known ground truth: how much accuracy is lost in a dividing population, and how much cell-cycle regression recovers (Fig. 5b) |
| The benefit of treatment holidays depends on a cost of resistance and on turnover (Fig. 8a-e) | Compartment models of adaptive therapy [@zhang2017adaptive; @strobl2021] | The same dependence from single cells whose resistance is a heritable expression state; cost and degree of protection become cell properties that lineage experiments measure; partial protection reverses the ranking even without a cost |
| Perturbing the retention of a resistant state changes the outcome (Fig. 4f) | Phenotypic-switching theory [@gunnarsson2020] and memory-disrupting compounds [@harmange2023] | The timing requirement (disruption must continue during exposure) and a lineage-resolved readout, the number of surviving clones |
| An intermediate dose can be optimal against drug-induced tolerance (Fig. 4h) | Reported for induced persisters [@corigliano2025] | Reproduced; the accompanying benefit of release periods does not appear here, which the text attributes to the decay of the induced state in this model rather than to a disagreement about data |
| MGMT expression selects glioblastoma cells under temozolomide (Fig. 8f-h) | Phenotypic selection with stable inheritance [@lasri2020] | Pharmacokinetics, the RTOG 0525 regimens, consumption of MGMT by the drug, and the fractionation that minimises the final population at fixed cumulative dose |
| Dose-dense temozolomide does not improve survival | Clinical result [@gilbert2013] | A mechanistic account of why, and identification of the rate at which the drug consumes MGMT as the measurement that separates the regimens |

# Supplementary Note 11: Identifiability and numerical robustness of the calibration

Six parameters were fitted to six measured fractions, so the fit cannot determine all six. Three checks bound what the calibration does and does not establish.

![](figS3_identifiability.png)

**Supplementary Fig. 3 | What the fate fractions determine.** **a**, Root-mean-square error of the simulated fate fractions at the training concentrations when one parameter is moved away from its fitted value over its full range and the others are held at theirs, plotted against the ratio to the fitted value; the dotted line marks one binomial standard error of the measurement above the best fit. **b**, The same error over the memory of the resistant state and the fraction of cells in it, with the death parameters held at their fitted values; red outlines mark the combinations that lie within one standard error of the best. **c**, Simulated fates at 13 µM with the integration step halved (mean ± s.d. of three seeds).

The death parameters are the better determined: the maximal death rate stays within one standard error of the best fit over about a twofold range, and the EC$_{50}$, Hill coefficient and growth-arrest IC$_{50}$ over narrower ranges still (Supplementary Fig. 3a). The two promoter switching rates are not determined: each can move by roughly sixfold with no penalty the data can detect. The two-dimensional slice shows why (Supplementary Fig. 3b). A short memory with many resistant cells and a long memory with few produce nearly the same fate fractions, so the error surface is flat along that direction, and a majority of the combinations tested, spanning one to nine generations of memory and 1% to 25% of cells resistant, fit within the noise of the measurement. The memory quoted in the main text is therefore set by the prior taken from the lineage correlations, not by the fate fractions, and the resistant fraction that accompanies it should not be read as a measurement.

Halving the integration step from 0.5 h to 0.25 h changes each simulated fate fraction by less than the seed-to-seed spread (Supplementary Fig. 3c), so the discretisation of the death hazard is not a source of error at this step size.

A profile in which the remaining parameters are re-optimised at each fixed memory would be preferable to a slice, but at seven seconds per objective evaluation the inner search could not be given enough budget to be reliable: it returned fits three times worse than the known solution at short memory, which reflects the optimiser rather than the data. The slice avoids that failure at the cost of conditioning on the fitted death parameters.
