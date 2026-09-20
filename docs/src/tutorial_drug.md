# Drug treatment and persisters

A [`Perturbation`](@ref) combines a dose schedule with the effects it drives.

```julia
model = telegraph_model(k_on = 0.005, k_off = 0.005, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
x0    = initial_state(model; G_off = 1)

pert = Perturbation(PulsedDose(1.0; on = 40.0, off = 20.0, start = 50.0);
    effects = [
        DeathHazard(h_max = 0.5, EC50 = 0.3, m = 2.0, protect = :protein, K = 60.0, q = 4.0),
        GrowthInhibition(IC50 = 2.0),
        GrowthCost(:protein; K = 60.0, q = 4.0, max_cost = 0.3),   # resistant cells grow slower
        RateModulation(:k_on, d -> 1 + 2d),     # drug-induced activation of the resistance gene
    ])
st  = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2)/20), control = FreeGrowth(max_cells = 20_000))
res = simulate_population(model, x0, 1000, (0.0, 200.0); settings = st, perturbation = pert, rng = Xoshiro(1))
res.popsize          # kill curve
res.dose             # dose applied at each record
```

`DeathHazard` implements state-dependent survival: the hazard is a Hill
function of the dose, reduced by the concentration of a protective species (a
resistance protein). Cells with a high resistance state before the drug
survive preferentially, which is the pre-existing-state mechanism of
non-genetic drug tolerance.

Two effects tie the drug to the cell cycle and to the protein that resists it.
[`CycleSensitivity`](@ref) makes the death hazard depend on where a cell is in
its division cycle, for agents whose lesions are converted into death during
replication:

```julia
CycleSensitivity(baseline = 0.2, center = 0.5, width = 0.15)
```

scales every `DeathHazard` by `baseline + (1 - baseline) exp(-((φ - center)/width)^2/2)`,
so killing peaks at cycle progress `center` and falls to a fraction `baseline`
of the peak elsewhere; `baseline = 1` recovers a cycle-independent hazard.

[`SuicideConsumption`](@ref) removes the protective protein stoichiometrically,
one molecule per lesion repaired, as happens to MGMT under an alkylating agent:

```julia
SuicideConsumption(:protein; k = 300.0, K_m = 150.0)
```

Because lesions form in proportion to the dose, the pool follows the
*cumulative* exposure rather than the peak concentration, which is what
distinguishes a fractionated schedule from a bolus of the same total dose.
Resynthesis is whatever the reaction model provides.

Schedules: [`ConstantDose`](@ref), [`PulsedDose`](@ref) (drug holidays),
[`PiecewiseDose`](@ref), [`BolusPK`](@ref) (one-compartment pharmacokinetics)
or any function wrapped in `FunctionDose`.

Genetic perturbations for causal benchmarks:

```julia
Perturbation(; gene_perturbations = [GenePerturbation(:k_tx, 0.1; fraction = 0.5, t_start = 20.0)])
```
knocks the parameter down in a random half of the cells (flag inherited by
daughters and recorded per cell).

Clonal statistics under drug:

```julia
ft = fluctuation_test(model, x0, sn -> mean(sn.counts[:, 4] ./ sn.volume .> 60); n_clones = 50, generations = 6, settings = st)
ft.ratio            # variance ratio ≫ 1 indicates heritable resistance states
```
