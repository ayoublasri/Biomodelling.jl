# Growing and dividing populations

[`simulate_population`](@ref) runs every cell of a population through the
stochastic kinetics of a model while the cell grows, replicates its genes,
divides and, optionally, dies. All options are collected in
[`PopulationSettings`](@ref):

```julia
st = PopulationSettings(
    dt            = 0.1,                         # global update / recording step
    kernel        = HybridSSATau(0.1),           # default: HybridSSATau(dt)
    growth        = ExponentialGrowth(log(2)/20; cv = 0.1),
    size_control  = Sizer(2.0; cv = 0.05),       # or Adder(1.0), AgeTimer(20.0)
    partitioning  = BinomialPartition(σ = 0.02), # or BetaBinomialPartition(ρ = 0.05)
    replication   = Replication(0.5),            # gene copies double mid-cycle
    control       = ConstantN(),                 # or FreeGrowth(), LogisticGrowth(K)
    record_every  = 1,
    track_lineage = true,
)
res = simulate_population(model, x0, 500, (0.0, 200.0); settings = st, rng = Xoshiro(1))
```

The result holds a record per `dt`: `res.counts[i]` (cells × species),
`res.volume[i]`, `res.ids[i]`, ages, generations, clone labels, gene copy
numbers, and `res.popsize`. Use [`snapshot`](@ref), [`final_snapshot`](@ref)
and [`concentrations`](@ref) to access them.

Division: molecules are partitioned binomially with the volume fraction
`f ~ Normal(0.5, σ)`; promoter groups are inherited by both daughters when
unreplicated and split one copy per daughter after replication; each daughter
draws its own division target and growth rate.

Population control: `ConstantN()` replaces a random cell by every second
daughter (the original `exponential_growth` behaviour), `FreeGrowth()` keeps
all cells (subsampling above `max_cells` while tracking the true size),
`LogisticGrowth(K)` adds a density-dependent death hazard.

Starting from a snapshot: pass an `N0 × nspecies` matrix of per-founder states as
`x0` (for example `final_snapshot(previous).counts`) together with `V0 = ...volume`.

Reproducibility: every cell carries its own random number generator seeded from
`rng`, so results do not depend on the number of threads.

## Lineage records

`res.lineage` is a [`LineageTable`](@ref) with one row per cell (`id`, `parent`,
`clone`, `generation`, birth and end times, fate, volume and state at birth and
at division). From it:

```julia
heritability(res, :protein; relation = :mother_daughter)   # (r, n)
lineage_autocorrelation(res, :protein; max_generations = 6)
memory_timescale(res, :protein)                            # in generations and time units
follow_lineage(res, res.ids[end][1])                       # a single-cell trace
newick(res.lineage; founder = 1, t_end = 200.0)
```

## Founders and promoter states

When `x0` is a single state vector, every founder starts from it but the promoter
states are drawn uniformly at random (`randomize_promoters = true`), which is the
stationary distribution only for symmetric switching. For slow, asymmetric
switching, initialise the founders explicitly (a founders × species matrix, or the
snapshot of a burn-in run) so that the resistant fraction starts at its
stationary value: the memory of a state with `k_on + k_off = 1/3000` per hour is
3000 hours, longer than most burn-ins.
