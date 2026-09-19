# Migration from v1 (v0.3)

The v1 API still works through a deprecated compatibility layer:

| v1 | v2 |
|---|---|
| `Donne(model, initiale, T, tau, NoC, growth)` | `ReactionModel(reactions; params)` + `initial_state` + `PopulationSettings` |
| NamedTuple reactions `(name, rate, reactants, products, coeff_rea, coeff_pro)` | [`Reaction`](@ref)`(name, reactants, products, kinetics)` |
| names containing `"act"`, `"inhib"`, `"comb_*"` | [`Hill`](@ref) kinetics with explicit regulators |
| species named `*on*` / `*off*` | promoter groups (`promoters = [[:G_off, :G_on]]`) |
| `ssa(data)` | `simulate(model, x0, tspan; kernel = DirectSSA())` |
| `tauleap`, `tauleapswitch`, `non_negative_Poisson_tauleap` | `TauLeap(τ)`, `HybridSSATau(τ)` |
| `adaptive_tauleap` | `AdaptiveTauLeap()` |
| `exponential_growth(data, div_noise, alg, Ni)` | `simulate_population(model, x0, N, tspan; settings)` |
| `trans_index` | `volume = :proportional` on the reaction |
| `random_network`, `network_generator` | [`random_grn`](@ref) |
| `comp_g`, `comp_tau`, `HO_reaction`, ... | internal to `AdaptiveTauLeap` |

Differences worth knowing:

* Species vectors no longer contain a `:NULL` entry (the compatibility layer
  adds a zero column to reproduce the old shapes).
* Output arrays are `Int`; volumes are `Float64`.
* The population loop records a [`LineageTable`](@ref) by default.
* Reactions may have any number of reactants and products.
