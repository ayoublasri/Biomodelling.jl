# Single cells

## Defining a model

A model is a list of [`Reaction`](@ref)s. Reactants and products are species
symbols (or `Symbol => coefficient` pairs); kinetics are [`MassAction`](@ref),
[`Hill`](@ref) or [`Custom`](@ref). Parameter values are given by name.

```julia
using Biomodelling
rx = [
    Reaction("transcription", [], [:mRNA], MassAction(:k_tx); volume = :proportional),
    Reaction("mRNA decay",    [:mRNA], [], MassAction(:γ_m)),
    Reaction("translation",   [:mRNA], [:mRNA, :P], MassAction(:k_tl); volume = :none),
    Reaction("protein decay", [:P], [], MassAction(:γ_p)),
]
model = ReactionModel(rx; params = (k_tx = 5.0, γ_m = 1.0, k_tl = 4.0, γ_p = 0.2))
```

The `volume` keyword controls how a propensity scales with the cell volume `V`:
`:auto` (default) gives `V^(1 - order)`, so zero-order production is
proportional to `V` and bimolecular reactions to `1/V`; `:proportional`,
`:inverse` and `:none` override it.

Regulation uses concentration-based Hill functions:

```julia
Reaction("regulated tx", [], [:mRNA],
         Hill(:k_tx; activators = [:A], inhibitors = [:R], K = [2.0, :K_R], n = 2.0, basal = 0.05))
```

Promoter states are ordinary species declared as a promoter group, so that
they are inherited (not partitioned) at division and doubled at replication:

```julia
tm = telegraph_model(k_on = 0.2, k_off = 0.5, k_tx = 20.0, k_dm = 1.0)   # G_off ⇄ G_on → mRNA
tm.promoter_groups                                                       # [[1, 2]]
```

## Simulating

```julia
x0 = initial_state(tm; G_off = 1)
tr = simulate(tm, x0, (0.0, 50.0); kernel = DirectSSA(), saveat = 0:0.5:50, rng = Xoshiro(1))
tr[:mRNA]                     # column of counts at the save times
E  = ensemble_final(tm, x0, 30.0, 10_000)   # final states of 10 000 independent cells
```

Kernels: `DirectSSA()` (exact), `TauLeap(τ)` (fixed step, errors on negative
populations), `HybridSSATau(τ)` (fixed step with exact fallback), and
`AdaptiveTauLeap(; ε, n_critical)`.

## Checking against theory

```julia
pmf = telegraph_pmf(0:60, 0.2, 0.5, 20.0, 1.0)     # Beta-Poisson stationary law
```
