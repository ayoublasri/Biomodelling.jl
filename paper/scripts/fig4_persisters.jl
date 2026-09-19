# Figure 4: drug response and persisters from heritable resistance states.
include(joinpath(@__DIR__, "common.jl"))

const λ = log(2) / 20
const T_DRUG = 60.0

"""Resistance gene R (telegraph promoter) → mRNA → protective protein P. Optional pretreatment window in which
promoter switching is accelerated `pre_factor`-fold (a memory-disrupting modulator)."""
function resistance_model(; k_on = 0.01, k_off = 0.01, pre_window = nothing, pre_factor = 20.0)
    if pre_window === nothing
        return telegraph_model(; k_on, k_off, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2, gene = :R, mrna = :mRNA, protein = :P)
    end
    t0, t1 = pre_window
    boost(t) = (t0 <= t < t1) ? pre_factor : 1.0
    rx = [Reaction("activation", [:R_off], [:R_on], Custom((x, V, θ, t) -> θ[1] * boost(t); params = [:k_on])),
          Reaction("inactivation", [:R_on], [:R_off], Custom((x, V, θ, t) -> θ[1] * boost(t); params = [:k_off])),
          Reaction("transcription", [:R_on], [:R_on, :mRNA], MassAction(:k_tx); volume = :proportional),
          Reaction("mRNA degradation", [:mRNA], [], MassAction(:k_dm)),
          Reaction("translation", [:mRNA], [:mRNA, :P], MassAction(:k_tl); volume = :none),
          Reaction("protein degradation", [:P], [], MassAction(:k_dp))]
    ReactionModel(rx; params = (k_on = k_on, k_off = k_off, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2), promoters = [[:R_off, :R_on]])
end
x0(m) = initial_state(m; R_off = 1)
death(; h_max = 0.5) = DeathHazard(h_max = h_max, EC50 = 0.5, m = 2.0, protect = :P, K = 150.0, q = 4.0)
settings(; N = 5000, kw...) = PopulationSettings(; dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05),
                                                  partitioning = BinomialPartition(σ = 0.02), control = FreeGrowth(max_cells = N), kw...)
pi(m) = speciesindex(m, :P)
high_fraction(sn, m; thr = 150.0) = mean(sn.counts[:, pi(m)] ./ sn.volume .> thr)

# (a) population trajectories under three schedules
mem = resistance_model()
rows = Any[]
for (tag, sched) in (("continuous_1.0", PiecewiseDose([0.0, T_DRUG], [0.0, 1.0])),
                     ("pulsed_1.0_on20_off10", PulsedDose(1.0; on = 20.0, off = 10.0, start = T_DRUG)),
                     ("continuous_0.4", PiecewiseDose([0.0, T_DRUG], [0.0, 0.4])))
    r = simulate_population(mem, x0(mem), 500, (0.0, 200.0); settings = settings(record_every = 5), perturbation = Perturbation(sched; effects = [death()]), rng = Xoshiro(1))
    for i in eachindex(r.t)
        push!(rows, [tag, r.t[i], r.popsize[i], r.dose[i], high_fraction(snapshot(r, i), mem)])
    end
    @printf("(a) %-24s final N = %.0f\n", tag, r.popsize[end])
end
save_csv("fig4a_trajectories.csv", ["schedule", "t", "popsize", "dose", "high_fraction"], permutedims(reduce(hcat, rows)))

# (b)+(c) kill curves, decay rates and single-cell division / death times versus dose
rows = Any[]; rows_c = Any[]; rows_ct = Any[]
for d in (0.0, 0.25, 0.5, 1.0, 2.0)
    r = simulate_population(mem, x0(mem), 500, (0.0, 140.0); settings = settings(record_every = 5), perturbation = Perturbation(PiecewiseDose([0.0, T_DRUG], [0.0, d]); effects = [death()]), rng = Xoshiro(2))
    i0 = argmin(abs.(r.t .- T_DRUG))
    for i in i0:length(r.t)
        push!(rows, [d, r.t[i] - T_DRUG, r.popsize[i] / r.popsize[i0]])
    end
    i30 = argmin(abs.(r.t .- (T_DRUG + 30)))
    decay = -(log(r.popsize[i30]) - log(r.popsize[i0])) / 30
    lt = r.lineage
    div = [lt.end_time[i] - lt.birth_time[i] for i in eachindex(lt.id) if lt.fate[i] == :divided && lt.birth_time[i] >= T_DRUG]
    dead = [lt.end_time[i] - max(lt.birth_time[i], T_DRUG) for i in eachindex(lt.id) if lt.fate[i] == :died]
    push!(rows_c, [d, decay, length(div), isempty(div) ? NaN : mean(div), isempty(div) ? NaN : std(div) / mean(div),
                   length(dead), isempty(dead) ? NaN : mean(dead), isempty(dead) ? NaN : std(dead) / mean(dead)])
    for v in div; push!(rows_ct, [d, "division", v]); end
    for v in dead; push!(rows_ct, [d, "death", v]); end
    @printf("(c) dose %.2f decay rate %.4f  division time %.2f (n=%d)  death time %.2f (n=%d)\n", d, decay, mean(div), length(div), isempty(dead) ? NaN : mean(dead), length(dead))
end
save_csv("fig4b_killcurves.csv", ["dose", "t_since_drug", "surviving_fraction"], permutedims(reduce(hcat, rows)))
save_csv("fig4c_decay_vs_dose.csv", ["dose", "decay_rate", "n_divisions", "division_time_mean", "division_time_cv", "n_deaths", "death_time_mean", "death_time_cv"], permutedims(reduce(hcat, rows_c)))
save_csv("fig4c_times.csv", ["dose", "event", "time"], permutedims(reduce(hcat, rows_ct)))

# (d) fate correlations between related cells (memory gene versus fast-switching control)
function fate_concordance(lt, pairs, t_from, t_to)
    same = 0; n = 0; deaths = 0
    for (a, b) in pairs
        (t_from <= lt.birth_time[a] < t_to && t_from <= lt.birth_time[b] < t_to) || continue
        fa = lt.fate[a] == :died; fb = lt.fate[b] == :died
        n += 1; same += (fa == fb); deaths += fa + fb
    end
    p = n == 0 ? NaN : deaths / (2n)
    (concordance = n == 0 ? NaN : same / n, expected = p^2 + (1 - p)^2, n = n, death_fraction = p)
end
rows = Any[]
for (tag, m) in (("memory", resistance_model(k_on = 0.01, k_off = 0.01)), ("fast", resistance_model(k_on = 1.0, k_off = 1.0)))
    r = simulate_population(m, x0(m), 800, (0.0, 110.0); settings = settings(N = 20000, record_every = 100), perturbation = Perturbation(PiecewiseDose([0.0, T_DRUG], [0.0, 1.0]); effects = [death()]), rng = Xoshiro(3))
    lt = r.lineage
    for (rel, pairs) in (("sisters", sister_pairs(lt)), ("cousins", Biomodelling.cousin_pairs(lt)))
        fc = fate_concordance(lt, pairs, T_DRUG - 25.0, T_DRUG)
        push!(rows, [tag, rel, fc.concordance, fc.expected, fc.n, fc.death_fraction])
        @printf("(d) %-6s %-8s concordance %.3f expected %.3f (n=%d)\n", tag, rel, fc.concordance, fc.expected, fc.n)
    end
end
save_csv("fig4d_fate_concordance.csv", ["model", "relation", "concordance", "expected_independent", "n_pairs", "death_fraction"], permutedims(reduce(hcat, rows)))

# (e) clone (barcode) diversity before and after drug: pre-existing versus drug-induced tolerance
effective_clones(clones) = (p = values(countmap(clones)) ./ length(clones); exp(-sum(p .* log.(p))))
rows = Any[]
induced = resistance_model(k_on = 0.0005, k_off = 0.01)
for (tag, m, effects) in (("pre_existing", mem, [death()]), ("drug_induced", induced, [death(), RateModulation(:k_on, d -> 1 + 200d)]))
    for seed in 1:3
        r = simulate_population(m, x0(m), 500, (0.0, 120.0); settings = settings(N = 20000, record_every = 10), perturbation = Perturbation(PiecewiseDose([0.0, T_DRUG], [0.0, 1.0]); effects = effects), rng = Xoshiro(10 + seed))
        pre = snapshot(r; t = T_DRUG); post = final_snapshot(r)
        push!(rows, [tag, seed, effective_clones(pre.clone), effective_clones(post.clone), size(pre.counts, 1), size(post.counts, 1), high_fraction(pre, m), high_fraction(post, m)])
    end
end
save_csv("fig4e_clone_diversity.csv", ["model", "seed", "effective_clones_before", "effective_clones_after", "cells_before", "cells_after", "high_fraction_before", "high_fraction_after"], permutedims(reduce(hcat, rows)))

# (f) schedule optimisation: release period × dose, pre-existing versus drug-induced tolerance
rows = Any[]
for (tag, m, extra) in (("pre_existing", mem, DrugEffect[]), ("drug_induced", induced, [RateModulation(:k_on, d -> 1 + 200d)]))
    for off in (0.0, 5.0, 10.0, 20.0, 40.0), d in (0.25, 0.5, 1.0, 2.0)
        sched = off == 0 ? PiecewiseDose([0.0, T_DRUG], [0.0, d]) : PulsedDose(d; on = 20.0, off = off, start = T_DRUG)
        r = simulate_population(m, x0(m), 300, (0.0, 260.0); settings = settings(N = 6000, record_every = 50), perturbation = Perturbation(sched; effects = vcat([death()], extra)), rng = Xoshiro(4))
        i0 = argmin(abs.(r.t .- T_DRUG))
        fitness = (log(max(r.popsize[end], 1e-9)) - log(r.popsize[i0])) / (260 - T_DRUG)
        exposure = sum(r.dose[i0:end]) / length(r.dose[i0:end])
        push!(rows, [tag, off, d, fitness, exposure, r.popsize[end]])
        @printf("(f) %-12s off %4.0f dose %.2f  fitness %.4f  final N %.0f\n", tag, off, d, fitness, r.popsize[end])
    end
end
save_csv("fig4f_schedules.csv", ["model", "release_period", "dose", "long_term_growth_rate", "mean_exposure", "final_popsize"], permutedims(reduce(hcat, rows)))

# (g) memory disruption: accelerate switching before the drug (pretreatment) and count surviving clones
rows = Any[]
for seed in 1:5
    for (tag, m) in (("no_pretreatment", resistance_model()), ("pretreatment", resistance_model(pre_window = (20.0, 60.0), pre_factor = 20.0)))
        r = simulate_population(m, x0(m), 500, (0.0, 140.0); settings = settings(N = 20000, record_every = 100), perturbation = Perturbation(PiecewiseDose([0.0, T_DRUG], [0.0, 1.0]); effects = [death()]), rng = Xoshiro(20 + seed))
        pre = snapshot(r; t = T_DRUG); post = final_snapshot(r)
        push!(rows, [tag, seed, length(unique(post.clone)), size(post.counts, 1), high_fraction(pre, m), high_fraction(post, m)])
        @printf("(g) %-16s seed %d surviving clones %d cells %d (high fraction before drug %.3f)\n", tag, seed, length(unique(post.clone)), size(post.counts, 1), high_fraction(pre, m))
    end
end
save_csv("fig4g_memory_disruption.csv", ["treatment", "seed", "surviving_clones", "surviving_cells", "high_fraction_before_drug", "high_fraction_after"], permutedims(reduce(hcat, rows)))

# (h) MGMT-like phenotypic selection: transient drug pulse enriches high expressers; enrichment persists with slow switching
rows = Any[]
for (tag, m) in (("slow_switching", resistance_model(k_on = 0.005, k_off = 0.005)), ("fast_switching", resistance_model(k_on = 0.2, k_off = 0.2)))
    r = simulate_population(m, x0(m), 500, (0.0, 260.0); settings = settings(N = 4000, record_every = 10), perturbation = Perturbation(PiecewiseDose([0.0, 60.0, 100.0], [0.0, 1.0, 0.0]); effects = [death()]), rng = Xoshiro(5))
    for i in eachindex(r.t)
        sn = snapshot(r, i)
        push!(rows, [tag, r.t[i], r.dose[i], r.popsize[i], mean(sn.counts[:, pi(m)] ./ sn.volume), high_fraction(sn, m)])
    end
end
save_csv("fig4h_mgmt.csv", ["model", "t", "dose", "popsize", "mean_P_concentration", "high_fraction"], permutedims(reduce(hcat, rows)))
println("fig4 done")
