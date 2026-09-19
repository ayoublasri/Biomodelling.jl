# Dose and schedule optimisation

The perturbation layer turns any population model into a testbed for treatment
schedules. This tutorial calibrates a persister model to laboratory-style
targets, compares clinical schedules, and searches for better ones.

## Clinical schedules

Schedules are built from a few primitives. Time units are the units of your
model's rate constants (hours below).

```julia
using Biomodelling, Random

# five days on, 23 off, six cycles (a 5/28 regimen), constant exposure during the days on
standard = PulsedDose(1.0; on = 5 * 24, off = 23 * 24, start = 0.0, cycles = 6)

# daily oral boluses with a 2.1 h elimination half-life on days 1-21 of 28-day cycles
dense = daily_boluses(0.5, cycle_days(21, 28, 6), log(2) / 2.1)

# adaptive therapy: stop when the tumour has shrunk to 50 % of its initial size,
# restart when it regrows to its initial size (Zhang et al. 2017)
adaptive = AdaptiveDose(1.0; on_above = 1.0, off_below = 0.5)
```

`AdaptiveDose` is evaluated with the population size, `dose(schedule, t, N)`;
`simulate_population` does this automatically.

## Outcomes

```julia
model = telegraph_model(k_on = 0.01, k_off = 0.05, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2,
                        gene = :R, mrna = :mRNA, protein = :P)
x0 = initial_state(model; R_off = 1)
st = PopulationSettings(dt = 0.5, growth = ExponentialGrowth(log(2) / 24), size_control = Sizer(2.0; cv = 0.1),
                        control = FreeGrowth(max_cells = 5000), record_every = 4)
death = DeathHazard(h_max = 0.05, EC50 = 0.5, m = 2.0, protect = :P, K = 150.0, q = 4.0)
run(sched; seed = 1) = simulate_population(model, x0, 300, (0.0, 24 * 60.0); settings = st,
                                           perturbation = Perturbation(sched; effects = [death]), rng = Xoshiro(seed))

r = run(standard)
net_growth_rate(r)                     # per hour, over the whole course
log_kill(r)                            # depth of response in log10 units
time_to_progression(r; threshold = 1.2) # RECIST-like progression from the nadir
cumulative_dose(r)                     # dose actually applied (also for adaptive schedules)
```

## Optimising a schedule under a dose constraint

`optimize_schedule` minimises an objective over a box of schedule parameters.
The objective receives the parameters and a seed; passing the same seeds to
every candidate compares schedules with common random numbers.

```julia
budget = cumulative_dose(standard, 0.0, 24 * 168.0)      # keep the cumulative dose of the standard regimen

function objective(θ, seed)
    period, duty, d = θ
    sched = PulsedDose(d; on = duty * period, off = (1 - duty) * period)
    over = max(0.0, cumulative_dose(sched, 0.0, 24 * 168.0) / budget - 1)
    -time_to_progression(run(sched; seed = seed))[1] / 24 + 1e3 * over    # days, with a penalty
end

opt = optimize_schedule(objective, [24.0, 0.1, 0.25], [24 * 28.0, 1.0, 2.0];
                        names = [:period_h, :duty, :dose], n_grid = 4, seeds = 1:2, maxiter = 30)
opt.params, opt.value
```

The result also carries `opt.table`, every candidate evaluated, which is what
Figure 8 of the paper plots. Objectives are cheap to change: extinction
probability over replicates (`extinction_probability`), final population size,
or a weighted sum with toxicity proxies such as the peak or cumulative dose.

## Calibrating to laboratory data first

Schedule predictions are only as good as the model behind them. The paper's
workflow (Figure 7) fits the persister model to time-lapse fate fractions at
one drug concentration with `abc_smc` and then checks the fitted model against
held-out concentrations, the biphasic kill curve and lineage correlations
before any schedule is optimised; `paper/scripts/fig7_calibration.jl` is a
template for that protocol.
