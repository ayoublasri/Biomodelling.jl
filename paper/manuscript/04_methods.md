# Methods

## Reaction model and propensities

A model is a set of species and reaction channels. Species are molecules (partitioned at division) or promoter states (inherited). For a reaction with reactant counts $x_i$, stoichiometric coefficients $\nu_i$, rate constant $k$ and volume exponent $e$, the mass-action propensity is
$$a = k\,V^{e}\prod_i \binom{x_i}{\nu_i},$$
where by default $e = 1-\sum_i \nu_i$, so that zero-order production scales with volume, first-order reactions are volume independent and bimolecular reactions scale as $1/V$; this keeps concentrations invariant under growth. Regulated reactions multiply this by $b + (1-b)F(\mathbf{c})$, where $b$ is a basal fraction and $F$ combines Hill terms of the regulator concentrations $c_j = x_j/V$, $h_j = c_j^{n_j}/(K_j^{n_j}+c_j^{n_j})$ for activators and $1-h_j$ for inhibitors, either as a product (AND logic) or as $1-\prod_j(1-h_j)$ (OR logic). Promoter switching is represented by explicit promoter species (for instance $G_{\mathrm{off}} \rightleftharpoons G_{\mathrm{on}}$) grouped into promoter groups. Custom propensities can be supplied as functions of counts, volume, parameters and time.

## Stochastic simulation kernels

Within one update step the cell volume and gene copy number are held fixed and the reaction network is advanced by one of four kernels: the direct method [@gillespie1977] with propensities updated through a species-reaction dependency graph; fixed-step Poisson tau-leaping [@gillespie2001]; a hybrid scheme in which a leap that would drive any species negative is rejected and the interval is instead simulated exactly (this is the exact version of the "switching" scheme of version 1); and adaptive tau-leaping with the step-size selection of @cao2006 and the critical-reaction partition of @cao2005, falling back to a fixed number of exact events when the proposed leap is shorter than $10/a_0$. Every kernel takes an explicit random number generator.

## Cells, growth, replication and division

Each cell carries its volume $V$, age, generation, lineage identifiers, gene copy state, growth rate, molecule counts, a cell-specific parameter vector and its own random number generator (seeded from the parent's generator, which makes results independent of the number of threads). Volume grows exponentially, $\dot V=\lambda V$, with an optional log-normal spread of $\lambda$ between cells. Division is triggered by a sizer ($V \ge V_{\mathrm{div}}$), an adder ($V - V_{\mathrm{birth}} \ge \Delta$) or a timer (age $\ge T$), each with a per-cell noise factor drawn at birth. At division the first daughter inherits a volume fraction $f\sim\mathcal N(0.5,\sigma^2)$ (clipped to $[0.05,0.95]$) and each molecule is assigned to it independently with probability $f$ (binomial partitioning); a beta-binomial option models over-dispersed partitioning. Promoter groups are inherited by both daughters when unreplicated and split one copy per daughter, in random assignment, when replicated. Gene replication, if enabled, doubles the promoter counts (and the propensity of reactions flagged as copy-number dependent) once the cell has completed a set fraction of its cycle, and is reset at division. Population control is either constant size (every second daughter replaces a uniformly chosen cell, the behaviour of version 1), free growth (all daughters kept, with uniform subsampling above a cap and book-keeping of the true size) or logistic (free growth with an additional hazard $\lambda N/K$).

## Drug and perturbation layer

A perturbation combines a dose schedule $d(t)$ (constant, pulsed, piecewise, one-compartment bolus pharmacokinetics or an arbitrary function) with effects evaluated for every cell at every step: a death hazard
$$h = h_{\max}\frac{d^m}{\mathrm{EC}_{50}^m + d^m}\cdot\frac{K^q}{K^q + c_P^{\,q}},$$
where $c_P$ is the concentration of a protective species (survival probability over a step is $e^{-h\,\Delta t}$); growth inhibition by $1/(1+(d/\mathrm{IC}_{50})^m)$; multiplicative modulation of any parameter by a function of the dose; and genetic perturbations (knockdown, knockout, over-expression) of any parameter in a random fraction of cells, the flag being inherited by daughters.

## Lineage statistics

The lineage table records, for every cell, its parent, founder (clone), generation, birth and end times, fate (divided, died, removed, alive) and its volume and state at birth and at the end of its life. From it we compute Pearson correlations of expression at division between mothers and daughters, sisters and cousins; the lineage autocorrelation $r_g$ between a cell and its ancestor $g$ generations earlier, and the memory timescale $\tau$ from $r_g = e^{-g/\tau}$; single-lineage traces by following the ancestral line of a surviving cell; the coefficient of variation squared of concentration across the population versus along single lineages; Luria-Delbrück fluctuation tests (variance of a phenotype across clones grown from single cells compared with the binomial expectation); and MemorySeq-style clonal scores [@shaffer2020], the variance of clone-mean expression divided by its permutation-null expectation.

## Observation models

True counts are converted to observations by a scRNA-seq model in which each cell captures every molecule independently with a Beta-distributed efficiency (mean and coefficient of variation set by the protocol), optional multinomial down-sampling to a Poisson-distributed depth, optional additional dropout and batch scaling; an smFISH model with a fixed detection efficiency; and a time-lapse model that samples a species along lineages at a fixed interval with additive or multiplicative measurement noise. Because counts scale with volume, library size correlates with cell size, reproducing the transcriptome-size confound of real data. Outputs are written as CSV or as AnnData-compatible HDF5 files.

## Inference

Snapshot counts of a gene in non-dividing cells follow the Beta-Poisson law of the telegraph model [@peccoud1995], which we evaluate through Kummer's confluent hypergeometric function with rescaled series summation, and maximise by a grid-initialised Nelder-Mead search in log-parameter space. For any model, including populations, we implement sequential Monte Carlo approximate Bayesian computation [@toni2009] with an adaptive tolerance schedule (the $\alpha$ quantile of the previous generation's distances), Gaussian perturbation kernels of variance twice the weighted population variance, importance weights against the prior, and parallel evaluation of the user-supplied distance function.

## Random regulatory networks

Random networks are built with a chosen number of activating and inhibiting edges under Erdős-Rényi or preferential-attachment (scale-free) sampling; regulation acts on the transcription rate of constitutive genes or on the promoter activation rate of telegraph genes through AND-combined Hill functions; per-gene rate constants can be drawn log-uniformly from ranges.

## Validation protocols

Exactness was tested by Kolmogorov-Smirnov distances between simulated stationary distributions and the analytic laws (Poisson for birth-death, Beta-Poisson for the telegraph model, negative binomial for bursty protein production [@friedman2006; @shahrezaei2008]), by agreement between kernels, and by a two-sample comparison against JumpProcesses.jl [@sciml] on the same model. Population-level properties were compared with the analytic expectations for cell-size scaling and concentration homeostasis [@bertaux2018], partitioning noise [@huh2011], lineage versus population statistics [@zhang2025] and cell-cycle-dependent burst parameters [@sukys2025]. Persister simulations were compared with the quantitative signatures reported by @iyer2025, @corigliano2025 and @harmange2023, and with the phenotypic-selection model of @lasri2020. All scripts, seeds and outputs are in the `paper/` directory of the repository.

## Software

Biomodelling.jl 2.0 is written in Julia (≥ 1.10) and depends on Distributions.jl, PoissonRandom.jl, QuadGK.jl, Optim.jl, SpecialFunctions.jl and StatsBase.jl; HDF5.jl is an optional extension. Benchmarks used NetworkInference.jl (PIDC, @chan2017), scikit-learn random forests re-implementing GENIE3 [@huynhthu2010], MAGIC [@vandijk2018] and kNN-smoothing [@wagner2018]. The test suite (Aqua quality checks and statistical tests) runs under GitHub Actions on Julia 1.10 and the latest release.
