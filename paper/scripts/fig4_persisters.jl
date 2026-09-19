# Figure 4: drug response and persisters from a rare, heritable resistance state.
# Usage: julia fig4_persisters.jl [a] [b] [d] [e] [f] [g] [h]   (default: all panels; b also produces c)
include(joinpath(@__DIR__, "common.jl"))
const parts = isempty(ARGS) ? ["a", "b", "d", "e", "f", "g", "h"] : ARGS

const λ = log(2) / 20
const T_DRUG = 40.0          # drug start in the treatment stage (after a burn-in of the founder population)

"""Resistance gene R (telegraph promoter) → mRNA → protective protein P. Optional pretreatment window in which
promoter switching is accelerated `pre_factor`-fold (a memory-disrupting modulator)."""
function resistance_model(; k_on = 0.002, k_off = 0.02, pre_window = nothing, pre_factor = 20.0)
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
death(; h_max = 0.15) = DeathHazard(h_max = h_max, EC50 = 0.5, m = 2.0, protect = :P, K = 150.0, q = 4.0)
growing(; N = 20000, kw...) = PopulationSettings(; dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05),
                                                  partitioning = BinomialPartition(σ = 0.02), control = FreeGrowth(max_cells = N), kw...)
const burn_settings = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05),
                                         partitioning = BinomialPartition(σ = 0.02), control = ConstantN(), track_lineage = false, record_every = 10_000)
"""Burn in a constant-size population so that promoter states, volumes and ages are at their stationary distribution."""
function burnin(m; N = 500, T = 300.0, seed = 0)
    r = simulate_population(m, x0(m), N, (0.0, T); settings = burn_settings, rng = Xoshiro(1000 + seed))
    sn = final_snapshot(r)
    (counts = sn.counts, volume = sn.volume)
end
"""Treatment-stage run started from a burnt-in founder population."""
function treat(m, pert; N = 500, T = 200.0, seed = 0, N_max = 20000, record_every = 5)
    b = burnin(m; N, seed)
    simulate_population(m, b.counts, N, (0.0, T); settings = growing(N = N_max, record_every = record_every), V0 = b.volume,
                        perturbation = pert, rng = Xoshiro(seed))
end
pidx(m) = speciesindex(m, :P)
high_fraction(sn, m; thr = 150.0) = mean(sn.counts[:, pidx(m)] ./ sn.volume .> thr)
from(t_on, d) = PiecewiseDose([0.0, t_on], [0.0, d])

mem = resistance_model()                                   # p_on = 0.002/0.022 ≈ 9%, memory time 1/(k_on+k_off) ≈ 45
fast = resistance_model(k_on = 0.05, k_off = 0.5)          # same p_on, memory time ≈ 1.8
induced = resistance_model(k_on = 0.0002, k_off = 0.02)   # p_on ≈ 1% before drug; the drug raises k_on 100-fold at dose 1
induction = RateModulation(:k_on, d -> 1 + 100d)

if "a" in parts
# (a) population trajectories under three schedules
rows = Any[]
for (tag, sched) in (("continuous_1.0", from(T_DRUG, 1.0)), ("pulsed_1.0_on20_off10", PulsedDose(1.0; on = 20.0, off = 10.0, start = T_DRUG)), ("continuous_0.4", from(T_DRUG, 0.4)))
    r = treat(mem, Perturbation(sched; effects = [death()]); seed = 1)
    for i in eachindex(r.t)
        push!(rows, [tag, r.t[i], r.popsize[i], r.dose[i], high_fraction(snapshot(r, i), mem)])
    end
    @printf("(a) %-24s N(drug) = %.0f  final N = %.0f\n", tag, r.popsize[argmin(abs.(r.t .- T_DRUG))], r.popsize[end])
end
save_csv("fig4a_trajectories.csv", ["schedule", "t", "popsize", "dose", "high_fraction"], permutedims(reduce(hcat, rows)))
end

if "b" in parts
# (b)+(c) kill curves, decay rates and single-cell division / death times versus dose
rows = Any[]; rows_c = Any[]; rows_ct = Any[]
for d in (0.0, 0.25, 0.5, 1.0, 2.0)
    r = treat(mem, Perturbation(from(T_DRUG, d); effects = [death()]); seed = 2, T = 140.0)
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
    @printf("(c) dose %.2f decay rate %.4f  division time %.2f (n=%d)  death time %.2f (n=%d)\n", d, decay, isempty(div) ? NaN : mean(div), length(div), isempty(dead) ? NaN : mean(dead), length(dead))
end
save_csv("fig4b_killcurves.csv", ["dose", "t_since_drug", "surviving_fraction"], permutedims(reduce(hcat, rows)))
save_csv("fig4c_decay_vs_dose.csv", ["dose", "decay_rate", "n_divisions", "division_time_mean", "division_time_cv", "n_deaths", "death_time_mean", "death_time_cv"], permutedims(reduce(hcat, rows_c)))
save_csv("fig4c_times.csv", ["dose", "event", "time"], permutedims(reduce(hcat, rows_ct)))
end

if "d" in parts
# (d) fate correlations between related cells (memory gene versus fast-switching control). A cell's fate is the
# fate of its lineage: it "survives" if any descendant is alive at the end of the run (colony formation), and
# "dies" if its whole subtree is extinct.
function lineage_survival(lt, alive_ids)
    cm = Biomodelling.children_map(lt)
    alive = Set(alive_ids)
    memo = Dict{Int,Bool}()
    function surv(id)
        haskey(memo, id) && return memo[id]
        v = id in alive || any(surv, get(cm, id, Int[]))
        memo[id] = v
    end
    surv
end
function fate_concordance(lt, pairs, t_from, t_to, surv)
    same = 0; n = 0; deaths = 0
    for (a, b) in pairs
        (t_from <= lt.birth_time[a] < t_to && t_from <= lt.birth_time[b] < t_to) || continue
        fa = !surv(a); fb = !surv(b)
        n += 1; same += (fa == fb); deaths += fa + fb
    end
    p = n == 0 ? NaN : deaths / (2n)
    (concordance = n == 0 ? NaN : same / n, expected = p^2 + (1 - p)^2, n = n, death_fraction = p)
end
rows = Any[]
for (tag, m) in (("memory", mem), ("fast", fast))
    r = treat(m, Perturbation(from(T_DRUG, 1.0); effects = [death()]); seed = 3, T = 100.0, N = 800, record_every = 100)
    lt = r.lineage
    surv = lineage_survival(lt, r.ids[end])
    for (rel, pairs) in (("sisters", sister_pairs(lt)), ("cousins", Biomodelling.cousin_pairs(lt)))
        fc = fate_concordance(lt, pairs, T_DRUG - 25.0, T_DRUG, surv)
        push!(rows, [tag, rel, fc.concordance, fc.expected, fc.n, fc.death_fraction])
        @printf("(d) %-6s %-8s concordance %.3f expected %.3f (n=%d, death fraction %.2f)\n", tag, rel, fc.concordance, fc.expected, fc.n, fc.death_fraction)
    end
end
save_csv("fig4d_fate_concordance.csv", ["model", "relation", "concordance", "expected_independent", "n_pairs", "death_fraction"], permutedims(reduce(hcat, rows)))
end

if "e" in parts
# (e) clone (barcode) diversity before and after drug: pre-existing versus drug-induced tolerance
effective_clones(clones) = (p = collect(values(countmap(clones))) ./ length(clones); exp(-sum(p .* log.(p))))
rows = Any[]
for (tag, m, effects) in (("pre_existing", mem, [death()]), ("drug_induced", induced, [death(), induction]))
    for seed in 1:3
        r = treat(m, Perturbation(from(T_DRUG, 1.0); effects = effects); seed = 10 + seed, T = 100.0, record_every = 10)
        pre = snapshot(r; t = T_DRUG); post = final_snapshot(r)
        push!(rows, [tag, seed, effective_clones(pre.clone), effective_clones(post.clone), size(pre.counts, 1), size(post.counts, 1), high_fraction(pre, m), high_fraction(post, m)])
        @printf("(e) %-12s seed %d clones %.0f → %.0f  cells %d → %d  high %.3f → %.3f\n", tag, seed, effective_clones(pre.clone), effective_clones(post.clone), size(pre.counts, 1), size(post.counts, 1), high_fraction(pre, m), high_fraction(post, m))
    end
end
save_csv("fig4e_clone_diversity.csv", ["model", "seed", "effective_clones_before", "effective_clones_after", "cells_before", "cells_after", "high_fraction_before", "high_fraction_after"], permutedims(reduce(hcat, rows)))
end

if "f" in parts
# (f) schedule optimisation: release period × dose for pre-existing tolerance (with and without a fitness cost) and drug-induced tolerance
cost = GrowthCost(:P; K = 150.0, q = 4.0, max_cost = 0.5)
rows = Any[]
for (tag, m, extra) in (("pre_existing", mem, DrugEffect[]), ("pre_existing_cost", mem, DrugEffect[cost]), ("drug_induced", induced, DrugEffect[induction]))
    for off in (0.0, 5.0, 10.0, 20.0, 40.0), d in (0.25, 0.5, 1.0, 2.0)
        sched = off == 0 ? from(T_DRUG, d) : PulsedDose(d; on = 20.0, off = off, start = T_DRUG)
        r = treat(m, Perturbation(sched; effects = vcat([death()], extra)); seed = 4, N = 300, T = 240.0, N_max = 6000, record_every = 50)
        i0 = argmin(abs.(r.t .- T_DRUG))
        fitness = (log(max(r.popsize[end], 1e-9)) - log(r.popsize[i0])) / (240 - T_DRUG)
        exposure = mean(r.dose[i0:end])
        push!(rows, [tag, off, d, fitness, exposure, r.popsize[end]])
        @printf("(f) %-18s off %4.0f dose %.2f  growth %.4f  final N %.0f\n", tag, off, d, fitness, r.popsize[end])
    end
end
save_csv("fig4f_schedules.csv", ["model", "release_period", "dose", "long_term_growth_rate", "mean_exposure", "final_popsize"], permutedims(reduce(hcat, rows)))
end

if "g" in parts
# (g) memory disruption: accelerate promoter switching 20-fold from two cell-cycle times before the drug, either
# stopping when the drug arrives (pretreatment only) or continuing throughout the exposure (co-treatment)
rows = Any[]
for seed in 1:5
    for (tag, m) in (("none", resistance_model()), ("before_drug", resistance_model(pre_window = (0.0, T_DRUG), pre_factor = 20.0)),
                     ("before_and_during", resistance_model(pre_window = (0.0, 120.0), pre_factor = 20.0)))
        r = treat(m, Perturbation(from(T_DRUG, 1.0); effects = [death()]); seed = 20 + seed, T = 120.0, record_every = 100)
        pre = snapshot(r; t = T_DRUG); post = final_snapshot(r)
        push!(rows, [tag, seed, length(unique(post.clone)), size(post.counts, 1), high_fraction(pre, m), high_fraction(post, m)])
        @printf("(g) %-18s seed %d surviving clones %d cells %d (high fraction before drug %.3f)\n", tag, seed, length(unique(post.clone)), size(post.counts, 1), high_fraction(pre, m))
    end
end
save_csv("fig4g_memory_disruption.csv", ["treatment", "seed", "surviving_clones", "surviving_cells", "high_fraction_before_drug", "high_fraction_after"], permutedims(reduce(hcat, rows)))
end

if "h" in parts
# (h) MGMT-like phenotypic selection: a drug pulse enriches high expressers; enrichment persists with slow switching
rows = Any[]
for (tag, m) in (("slow_switching", mem), ("fast_switching", fast))
    r = treat(m, Perturbation(PiecewiseDose([0.0, T_DRUG, T_DRUG + 40.0], [0.0, 1.0, 0.0]); effects = [death()]); seed = 5, T = 260.0, N_max = 4000, record_every = 10)
    for i in eachindex(r.t)
        sn = snapshot(r, i)
        push!(rows, [tag, r.t[i], r.dose[i], r.popsize[i], mean(sn.counts[:, pidx(m)] ./ sn.volume), high_fraction(sn, m)])
    end
end
save_csv("fig4h_mgmt.csv", ["model", "t", "dose", "popsize", "mean_P_concentration", "high_fraction"], permutedims(reduce(hcat, rows)))
end

println("fig4 done")
