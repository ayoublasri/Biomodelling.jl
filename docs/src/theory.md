# Theory notes

## Propensities

For a reaction with reactant counts `x_i`, coefficients `ν_i`, rate `k` and
volume exponent `e`, the mass-action propensity is
`k V^e ∏_i C(x_i, ν_i)`. With `volume = :auto`, `e = 1 − Σ ν_i`, which keeps
concentrations invariant under volume changes (production ∝ V, bimolecular
∝ 1/V). Regulated reactions multiply this by
`basal + (1 − basal) F(c)` with `F` the AND (product) or OR combination of Hill
terms of regulator concentrations `c = x/V`.

## Stationary distributions used for validation

* Birth-death: Poisson(k/γ).
* Telegraph (Peccoud and Ycart 1995): Beta-Poisson with `a = k_on/γ`,
  `b = k_off/γ`, `λ = k_tx/γ`; implemented in [`telegraph_pmf`](@ref) through
  Kummer's function.
* Bursty protein (Friedman, Cai and Xie 2006; Shahrezaei and Swain 2008):
  negative binomial with shape `a/γ` and success probability `1/(1 + b)`.

## Growth and division

Volume grows as `dV/dt = λV`; sizer, adder and timer rules trigger division;
molecules are partitioned binomially with fraction `f`; promoter groups are
inherited. Gene replication doubles promoter counts (and constitutive
transcription rates flagged `copy_number = true`) at a set cycle fraction. These
ingredients reproduce cell-size scaling of transcript numbers, concentration
homeostasis (Bertaux, Marguerat and Shahrezaei 2018), the partitioning-noise
contribution to variability (Huh and Paulsson 2011), and the difference between
lineage and population statistics (Zhang, Singh et al. 2025).

## Death under drug

The death hazard is `h_max d^m/(EC50^m + d^m) · K^q/(K^q + c^q)` with `c` the
concentration of the protective species; survival over a step `Δt` is
`exp(−h Δt)`.
