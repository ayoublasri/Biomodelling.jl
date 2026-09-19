---
title: "Supplementary information for: Biomodelling.jl 2.0"
author: "Ayoub Lasri and Marc Sturrock"
date: ""
---

# S1. Stationary laws used for validation

**Birth-death.** For $\emptyset \xrightarrow{k} X \xrightarrow{\gamma} \emptyset$ the stationary distribution is Poisson with mean $k/\gamma$.

**Telegraph model.** With activation $k_{\mathrm{on}}$, inactivation $k_{\mathrm{off}}$, transcription $k_{\mathrm{tx}}$ from the active state and degradation $\gamma$, and $a = k_{\mathrm{on}}/\gamma$, $b = k_{\mathrm{off}}/\gamma$, $\lambda = k_{\mathrm{tx}}/\gamma$, the stationary distribution is the Beta-Poisson mixture
$$P(n) = \frac{\lambda^n}{n!}\frac{\Gamma(a+n)\,\Gamma(a+b)}{\Gamma(a+b+n)\,\Gamma(a)}\,{}_1F_1(a+n;\,a+b+n;\,-\lambda),$$
whose mean is $\lambda a/(a+b)$ and whose variance is $\mu + \lambda^2 ab/[(a+b)^2(a+b+1)]$. We evaluate ${}_1F_1(\alpha;\beta;-\lambda)$ through Kummer's transformation $e^{-\lambda}{}_1F_1(\beta-\alpha;\beta;\lambda)$, whose series has positive terms; partial sums are rescaled by $10^{-200}$ whenever they exceed $10^{200}$ so that the log of the sum is accumulated without overflow.

**Bursty protein.** Protein produced in geometric bursts of mean size $b$ arriving at rate $a$ and degraded at rate $\gamma$ is negative binomial with shape $a/\gamma$ and success probability $1/(1+b)$. The package emulates bursts through an mRNA of lifetime $1/(200\gamma)$ translated at rate $200\gamma b$.

# S2. Kernel details and correctness tests

The direct method draws the waiting time from $\mathrm{Exp}(a_0)$ and the channel by linear search over the cumulative propensities; after each event only the propensities of reactions that depend on a changed species (the dependency graph) are recomputed, and the running total is refreshed every 512 events. Fixed-step tau-leaping draws $\mathrm{Poisson}(a_j\tau)$ firings for every channel; the strict variant errors when a species would become negative, the hybrid variant rejects the leap and simulates the interval with the direct method. Adaptive tau-leaping classifies as critical every channel within $n_c = 10$ firings of exhausting a reactant, selects $\tau$ from the non-critical channels through the bounds $\max(\varepsilon x_i/g_i, 1)/|\mu_i|$ and $\max(\varepsilon x_i/g_i, 1)^2/\sigma_i^2$, draws the time to the next critical event from $\mathrm{Exp}(a_0^{\mathrm{crit}})$, fires at most one critical channel, halves $\tau$ on rejection, and performs 100 exact events whenever $\tau < 10/a_0$.

Table S1 lists the tests of the package test suite that validate the kernels: Kolmogorov-Smirnov distances below 0.015 (direct) or 0.02 (approximate kernels) to the Poisson, Beta-Poisson and negative-binomial laws with 8 000 to 12 000 cells; agreement of the adaptive kernel with the direct method on a stiff dimerisation network within 5% in the mean; identical results with one and with several threads.

# S3. Parameters of all simulations

| Figure | Model | Parameters |
|---|---|---|
| 2a | birth-death | $k = 12$, $\gamma = 1$, $T = 15$ |
| 2b | telegraph | $(k_{\mathrm{on}}, k_{\mathrm{off}}, k_{\mathrm{tx}}) = (0.4, 0.6, 12)$ and $(0.05, 0.5, 40)$, $\gamma = 1$, $T = 40$ |
| 2c | bursty protein | $a = 1.5$, $b = 6$, $\gamma = 1$, $T = 25$ |
| 2d-e | telegraph | $(0.4, 0.6, 12)$, $n = 20\,000$ |
| 2f | random telegraph GRNs | $n_{\mathrm{act}} = n_{\mathrm{genes}}$, $n_{\mathrm{inh}} = n_{\mathrm{genes}}/2$, $\lambda = \ln 2/20$, $dt = 0.1$, $T = 20$ |
| 3 | telegraph gene with protein | $k_{\mathrm{tx}} = 30$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 4$, $k_{\mathrm{dp}} = 0.2$, $k_{\mathrm{on}} = k_{\mathrm{off}}$ as indicated; $\lambda = \ln 2/20$, sizer $V_{\mathrm{div}} = 2$ (cv 0.05), $\sigma_f = 0.02$, replication at 50% of the cycle where indicated |
| 3e | stable protein | $k_{\mathrm{on}} = k_{\mathrm{off}} = 5$, $k_{\mathrm{tx}} = 20$, $k_{\mathrm{dm}} = 1$, $k_{\mathrm{tl}} = 2$, $k_{\mathrm{dp}} = 0.05$ |
| 4 | resistance gene | as Figure 3 with $k_{\mathrm{on}} = k_{\mathrm{off}} = 0.01$ (memory) or 1 (fast); death $h_{\max} = 0.5$, $\mathrm{EC}_{50} = 0.5$, $m = 2$, $K = 150$, $q = 4$; drug from $t = 60$; drug-induced variant $k_{\mathrm{on}} = 0.0005$, $k_{\mathrm{off}} = 0.01$, $k_{\mathrm{on}} \to k_{\mathrm{on}}(1 + 200 d)$; pretreatment: switching $\times 20$ on $[20, 60)$ |
| 5a | 40 independent telegraph genes | $k_{\mathrm{on}} = k_{\mathrm{off}} \in [10^{-3}, 1]$ (log-spaced), $k_{\mathrm{tx}} = 20$, $k_{\mathrm{dm}} = 0.5$; 40 founders, 6 doublings |
| 5b-c | random GRN | 30 genes, 30 activations, 15 inhibitions, $k_{\mathrm{tx}} \in [5, 30]$, $K \in [2, 10]$, $n = 2$, basal 0.05; 1 500 cells; sequencing capture 0.15 (cv 0.3) |
| 5d | random GRN | 20 genes, 22 activations, 10 inhibitions; knockdown factor 0.05; 800 cells |
| 6a-b | telegraph | truth $(0.3, 0.6, 20)$; priors log-uniform on $[0.01, 10]$, $[0.01, 10]$, $[1, 200]$; 200 (a) / 120 (b) particles |
| 6c | resistance gene | truth $h_{\max} = 0.5$, $K = 150$; priors log-uniform on $[0.05, 5]$ and $[20, 1000]$; 100 particles |

# S4. Additional persister simulations

[Placeholder: dose-response of the death-time distributions (Fig. 4c data), drug-induced versus pre-existing tolerance under continuous dosing, sensitivity to the protection threshold $K$ and Hill coefficient $q$.]

# S5. Benchmark details

Network inference methods: absolute Pearson and Spearman correlations of $\log(1+x)$; PIDC as implemented in NetworkInference.jl; GENIE3 as random-forest importances (200 trees, `sqrt` features) fitted on $\log(1+x)$ with scikit-learn. Imputation: kNN-smoothing (one step, $k = 15$, 10 principal components of the Freeman-Tukey transformed, library-size normalised counts) and MAGIC with default parameters on square-root normalised counts. Metrics: area under the precision-recall and receiver-operating curves against the undirected (correlation, PIDC) or directed (GENIE3) ground-truth adjacency.

# S6. Inference diagnostics

[Placeholder: tolerance schedules and acceptance rates per generation for Figures 6a-c; posterior predictive checks of the summary statistics.]

# S7. Runtime

[Placeholder: Table of runtimes from Figure 2f and from `benchmarks/runtime.jl`; scaling is linear in the number of cells and approximately linear in the number of reactions.]

# S8. Migration from version 1

The v1 API (`Donne`, `ssa`, `tauleap`, `tauleapswitch`, `adaptive_tauleap`, `exponential_growth`, NamedTuple reactions) remains available through a deprecated compatibility layer that converts NamedTuple reactions to `ReactionModel`s (reactions whose names contain `act`, `inhib`, `comb_a` or `comb_i` become Hill kinetics with `rate = [k, n, K]`; species whose names contain `on`/`off` become promoter groups) and maps the algorithms onto the new kernels. Output shapes are preserved, including the leading `:NULL` column. See `docs/src/migration.md`.
