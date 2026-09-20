# Validation of the population and lineage layers against exact stationary solutions.
#
# Beentjes, Perez-Carrasco & Grima (Phys. Rev. E 2020, 101:032403) and Jia & Grima
# (iScience 2023, 26:105746) solve the stationary gene-expression distribution of a
# growing, dividing cell population exactly, and give different closed forms for a
# single lineage and for a population snapshot. This script generates the simulated
# counterparts of the analytically solvable members of that class; the exact solutions
# and the comparison are computed in `paper/scripts/exact_solutions.py`.
#
# Four models, each run in both modes:
#   constitutive  zero-order production, first-order decay, memoryless interdivision time
#   bursty        production in fixed bursts of `BURST` molecules, otherwise as above
#   telegraph     two-state promoter, transcription, decay, promoter inherited at division
#   replication   deterministic cycle, volume-scaled synthesis, gene replication at mid-cycle
#
# The first three have an exponentially distributed interdivision time, which makes the
# count process Markovian in the molecule number alone and yields a closed-form stationary
# law in each mode. The fourth has a deterministic cycle, for which the count distribution
# at every cycle phase is exactly Poisson and the two modes differ through the age density.
#
# Usage: julia fig9_validation.jl [quick]
include(joinpath(@__DIR__, "common.jl"))
const QUICK = "quick" in ARGS

const T_DIV = 1.0                       # mean interdivision time
const GAMMA = 1.0                       # first-order decay rate
const K_PROD = 20.0                     # zero-order production rate (molecules per unit time)
const BURST = 4                         # molecules per burst in the bursty model
const K_ON, K_OFF, K_TX = 0.5, 1.0, 40.0
const NMAX = 250                        # histogram range

# The step divides the interdivision time an exact (binary) number of times, so that the
# deterministic timer fires on the same step for every cell: with a step such as 0.002 the
# comparison `age >= T` falls on a floating-point boundary and half the cells divide one
# step late, which smears the cycle-phase lattice without changing the count distribution.
const DT = QUICK ? 1 / 64 : 1 / 512
const T_END = QUICK ? 6.0 : 12.0        # burn-in, in units of the mean interdivision time
const N_LINEAGE = QUICK ? 4_000 : 40_000
const N_POP_CAP = QUICK ? 4_000 : 15_000
const N_POP_REPS = QUICK ? 2 : 5
const DT_SCAN = QUICK ? [1 / 64, 1 / 128] : [1 / 64, 1 / 128, 1 / 256, 1 / 512]

# ---------------------------------------------------------------- models
# Every propensity is made volume-independent so that the comparison isolates division,
# partitioning and the lineage/snapshot weighting (the replication model below is the one
# that exercises volume scaling).
const M_CONST = ReactionModel([Reaction("production", [], [:M], MassAction(:k); volume = :none),
                               Reaction("decay", [:M], [], MassAction(:g); volume = :none)];
                              params = (k = K_PROD, g = GAMMA))
const M_BURST = ReactionModel([Reaction("burst", [], [:M => BURST], MassAction(:k); volume = :none),
                               Reaction("decay", [:M], [], MassAction(:g); volume = :none)];
                              params = (k = K_PROD / BURST, g = GAMMA))
const M_TELE = ReactionModel([Reaction("activation", [:G_off], [:G_on], MassAction(:k_on); volume = :none),
                              Reaction("inactivation", [:G_on], [:G_off], MassAction(:k_off); volume = :none),
                              Reaction("transcription", [:G_on], [:G_on, :M], MassAction(:k_tx); volume = :none),
                              Reaction("decay", [:M], [], MassAction(:g); volume = :none)];
                             params = (k_on = K_ON, k_off = K_OFF, k_tx = K_TX, g = GAMMA),
                             promoters = [[:G_off, :G_on]])
# Volume-scaled synthesis from a gene that replicates at mid-cycle, in a cell that grows
# exponentially and halves its volume at division (Jia & Grima 2023).
const LAM = log(2) / T_DIV
const REP_FRACTION = 0.5
const M_REPL = ReactionModel([Reaction("production", [], [:M], MassAction(:k); volume = :proportional, copy_number = true),
                              Reaction("decay", [:M], [], MassAction(:g); volume = :none)];
                             params = (k = K_PROD, g = GAMMA))

x0_of(m) = m === M_TELE ? initial_state(m; G_off = 1) : initial_state(m)

# ---------------------------------------------------------------- settings
"""Memoryless (exponentially distributed) interdivision time; volume plays no role."""
expo_settings(dt, control) = PopulationSettings(; dt = dt, kernel = DirectSSA(), growth = ExponentialGrowth(0.0),
                                                size_control = AgeTimer(T_DIV; cv = 1.0, dist = :gamma),
                                                partitioning = BinomialPartition(), control = control,
                                                track_lineage = false, record_every = 10_000_000, threads = true)
"""Deterministic cycle with exponential volume growth and gene replication at mid-cycle."""
repl_settings(dt, control) = PopulationSettings(; dt = dt, kernel = DirectSSA(), growth = ExponentialGrowth(LAM),
                                                size_control = AgeTimer(T_DIV), partitioning = BinomialPartition(),
                                                replication = Replication(REP_FRACTION), control = control,
                                                track_lineage = false, record_every = 10_000_000, threads = true)

hist_of(v) = (h = zeros(Int, NMAX + 1); for x in v; h[clamp(x, 0, NMAX)+1] += 1; end; h)

rows_counts = Any[]
push_hist!(case, mode, kernel, dt, rep, v) = for (i, c) in enumerate(hist_of(v))
    c == 0 && continue
    push!(rows_counts, [case, mode, kernel, dt, rep, i - 1, c])
end

"""Molecule counts of the final record of one run."""
function final_counts(m, settings, N0, T; seed, age0 = nothing, V0 = nothing)
    r = simulate_population(m, x0_of(m), N0, (0.0, T); settings = settings, rng = Xoshiro(seed),
                            age0 = age0, V0 = V0)
    sn = final_snapshot(r)
    (counts = sn.counts[:, speciesindex(m, :M)], age = sn.age, n = size(sn.counts, 1))
end

# ---------------------------------------------------------------- memoryless-timer models
for (case, m) in (("constitutive", M_CONST), ("bursty", M_BURST), ("telegraph", M_TELE))
    t0 = time()
    # single lineages: N_LINEAGE independent mother machines, so the sample is i.i.d.
    v = final_counts(m, expo_settings(DT, MotherMachine()), N_LINEAGE, T_END; seed = 11).counts
    push_hist!(case, "lineage", "DirectSSA", DT, 1, v)
    @printf("%-13s lineage    n=%6d  mean %.4f  var %.4f  (%.0f s)\n", case, length(v), mean(v), var(v), time() - t0)
    # population snapshots: replicate branching populations, uniformly subsampled at the cap
    for rep in 1:N_POP_REPS
        t1 = time()
        v = final_counts(m, expo_settings(DT, FreeGrowth(max_cells = N_POP_CAP)), 200, T_END; seed = 100rep + 3).counts
        push_hist!(case, "population", "DirectSSA", DT, rep, v)
        @printf("%-13s population n=%6d  mean %.4f  var %.4f  (rep %d, %.0f s)\n", case, length(v), mean(v), var(v), rep, time() - t1)
    end
end

# ---------------------------------------------------------------- step-size convergence
# The interdivision time is resolved to the update step `dt`, which is the only source of
# bias in the memoryless-timer comparison; this scan shows it vanishing.
for dt in DT_SCAN
    v = final_counts(M_CONST, expo_settings(dt, MotherMachine()), N_LINEAGE ÷ 2, T_END; seed = 21).counts
    push_hist!("constitutive_dt", "lineage", "DirectSSA", dt, 1, v)
    @printf("dt scan   dt=%.4f  n=%6d  mean %.4f\n", dt, length(v), mean(v))
end
# the hybrid tau-leaping kernel at the production step it is normally run with
let s = PopulationSettings(; dt = DT, kernel = HybridSSATau(DT), growth = ExponentialGrowth(0.0),
                           size_control = AgeTimer(T_DIV; cv = 1.0, dist = :gamma), partitioning = BinomialPartition(),
                           control = MotherMachine(), track_lineage = false, record_every = 10_000_000)
    v = final_counts(M_CONST, s, N_LINEAGE ÷ 2, T_END; seed = 31).counts
    push_hist!("constitutive", "lineage", "HybridSSATau", DT, 1, v)
    @printf("hybrid kernel lineage n=%6d  mean %.4f\n", length(v), mean(v))
end

# ---------------------------------------------------------------- deterministic cycle with replication
# Founder ages are placed uniformly on the update grid, which is the stationary age
# distribution of a single lineage; the population's age distribution is left to emerge
# from the branching dynamics and is recorded alongside the counts.
const J = round(Int, T_DIV / DT)
const T_REPL = QUICK ? 6.0 : 14.0                    # burn-in in generations
const T_REPL_HALF = T_REPL + 0.5                     # record half a cycle later, where the two age densities differ most
founder_ages(n) = [DT * ((i - 1) % J) for i in 1:n]
founder_vols(ages) = exp.(LAM .* ages)

rows_ages = Any[]
for (mode, control, n0, reps, tend) in (("lineage", MotherMachine(), QUICK ? 4_000 : 20_000, 1, T_REPL),
                                        ("population", FreeGrowth(max_cells = N_POP_CAP), 2J, N_POP_REPS, T_REPL_HALF))
    for rep in 1:reps
        t1 = time()
        a0 = founder_ages(n0)
        f = final_counts(M_REPL, repl_settings(DT, control), n0, tend; seed = 700rep + 5, age0 = a0, V0 = founder_vols(a0))
        push_hist!("replication", mode, "DirectSSA", DT, rep, f.counts)
        ah = zeros(Int, J)
        for a in f.age
            ah[clamp(round(Int, a / DT), 0, J - 1)+1] += 1
        end
        for (i, c) in enumerate(ah)
            c == 0 || push!(rows_ages, [mode, rep, i - 1, c])
        end
        @printf("replication   %-10s n=%6d  mean %.4f  var %.4f  (rep %d, %.0f s)\n", mode, f.n, mean(f.counts), var(f.counts), rep, time() - t1)
    end
end

save_csv("fig9_counts.csv", ["case", "mode", "kernel", "dt", "replicate", "n", "cells"], permutedims(reduce(hcat, rows_counts)))
save_csv("fig9_ages.csv", ["mode", "replicate", "age_step", "cells"], permutedims(reduce(hcat, rows_ages)))
save_kv("fig9_settings.csv", ["T_div" => T_DIV, "gamma" => GAMMA, "k_prod" => K_PROD, "burst" => BURST,
                              "k_on" => K_ON, "k_off" => K_OFF, "k_tx" => K_TX, "dt" => DT, "steps_per_cycle" => J,
                              "replication_fraction" => REP_FRACTION, "growth_rate" => LAM,
                              "t_end_expo" => T_END, "t_end_repl_lineage" => T_REPL,
                              "t_end_repl_population" => T_REPL_HALF, "n_lineage" => N_LINEAGE,
                              "pop_cap" => N_POP_CAP, "pop_replicates" => N_POP_REPS, "nmax" => NMAX])
println("fig9 done")
