# Figure 7: calibration of the persister model to laboratory time-lapse data (cisplatin, Iyer et al. 2025)
# and validation on held-out conditions. Time in hours, concentrations in µM.
# Usage: julia fig7_calibration.jl [calibrate] [validate]
include(joinpath(@__DIR__, "common.jl"))
const parts = isempty(ARGS) ? ["calibrate", "validate"] : ARGS
const DATA = joinpath(@__DIR__, "..", "data")

# ---------------------------------------------------------------- data (transcribed from the article; see paper/data/README.md)
raw, hdr = readdlm(joinpath(DATA, "iyer2025_u2os_fates.csv"), ','; header = true)
col(name) = findfirst(==(name), vec(hdr))
obs = Dict{Float64,NamedTuple}()
for i in 1:size(raw, 1)
    n = raw[i, col("cells_at_drug")]
    obs[Float64(raw[i, col("cisplatin_uM")])] = (died = raw[i, col("died")] / n, divided = raw[i, col("divided")] / n,
                                                  survived = raw[i, col("survived_without_dividing")] / n, n = n, label = String(raw[i, col("condition")]))
end
const TRAIN = [7.0, 13.0]          # Low and High are used for calibration
const HELDOUT = [10.0]             # Medium is predicted
const MEMORY_GENERATIONS = 2.5     # states inherited across 2-3 generations (HCT116 lineage correlations)

# ---------------------------------------------------------------- model (hours)
const CYCLE = 24.0                 # U2OS doubling time ≈ 24 h
const λ = log(2) / CYCLE
const T_FREE = 48.0                # two drug-free days of imaging
const T_DRUG = 72.0                # three days of cisplatin
model(k_on, k_off) = telegraph_model(; k_on, k_off, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2, gene = :R, mrna = :mRNA, protein = :P)
x0(m) = initial_state(m; R_off = 1)
effects(h_max, EC50, mh, IC50) = [DeathHazard(h_max = h_max, EC50 = EC50, m = mh, protect = :P, K = 150.0, q = 4.0),
                                  GrowthInhibition(IC50 = IC50, m = 2.0)]
burn_st = PopulationSettings(dt = 0.5, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.1), partitioning = BinomialPartition(σ = 0.02),
                             control = ConstantN(), track_lineage = false, record_every = 10_000)
treat_st(; record_every = 4) = PopulationSettings(dt = 0.5, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.1), partitioning = BinomialPartition(σ = 0.02),
                                                  control = FreeGrowth(max_cells = 20_000), track_lineage = true, record_every = record_every)
"""Burn in a constant-size population (promoter states, sizes and ages at stationarity) and treat it: `t_free` hours
without drug, then `t_drug` hours at concentration `d` (µM)."""
function experiment(θ, d; N = 300, seed = 1, t_free = T_FREE, t_drug = T_DRUG, record_every = 4)
    k_on, k_off, h_max, EC50, mh, IC50 = θ
    m = model(k_on, k_off)
    b = simulate_population(m, x0(m), N, (0.0, 10CYCLE); settings = burn_st, rng = Xoshiro(1000 + seed))
    sn = final_snapshot(b)
    pert = Perturbation(PiecewiseDose([0.0, t_free], [0.0, d]); effects = effects(h_max, EC50, mh, IC50))
    simulate_population(m, sn.counts, N, (0.0, t_free + t_drug); settings = treat_st(; record_every), V0 = sn.volume, perturbation = pert, rng = Xoshiro(seed))
end
"""Fates, over the drug window, of the cells present at drug addition: died, divided, or survived without dividing."""
function fates(r; t_on = T_FREE, t_end = T_FREE + T_DRUG)
    lt = r.lineage
    died = divided = survived = 0
    for i in eachindex(lt.id)
        lt.birth_time[i] <= t_on || continue
        (isnan(lt.end_time[i]) || lt.end_time[i] > t_on) || continue
        f = lt.fate[i]
        if f == :died && lt.end_time[i] <= t_end
            died += 1
        elseif f == :divided && lt.end_time[i] <= t_end
            divided += 1
        else
            survived += 1
        end
    end
    n = died + divided + survived
    (died = died / n, divided = divided / n, survived = survived / n, n = n)
end
fatevec(f) = [f.died, f.divided, f.survived]

# ---------------------------------------------------------------- calibration by simulated moments with common random numbers
names = [:k_on, :k_off, :h_max, :EC50, :m_h, :IC50]
lower = [1e-3, 1e-3, 0.005, 3.0, 1.0, 1.0]
upper = [0.3, 0.3, 0.5, 40.0, 6.0, 40.0]
function distance(θ, seed)
    d2 = 0.0
    for d in TRAIN
        f = fates(experiment(θ, d; seed))
        d2 += sum((fatevec(f) .- fatevec(obs[d])) .^ 2)
    end
    memgen = 1 / (θ[1] + θ[2]) / CYCLE
    sqrt(d2 / (3 * length(TRAIN))) + 0.05 * abs(log(memgen / MEMORY_GENERATIONS))
end
if "calibrate" in parts
    t0 = time()
    opt = optimize_schedule(distance, lower, upper; names, n_grid = 3, max_grid = 160, seeds = 1:1, maxiter = 80, rng = Xoshiro(7), verbose = false)
    @printf("calibration: %s  (%.0f s)\n", opt, time() - t0)
    save_kv("fig7_calibration.csv", vcat([string(n) => v for (n, v) in zip(names, opt.params)], ["distance" => opt.value, "n_evaluations" => length(opt.table),
            "memory_generations" => 1 / (opt.params[1] + opt.params[2]) / CYCLE, "p_on" => opt.params[1] / (opt.params[1] + opt.params[2])]))
    save_csv("fig7_calibration_table.csv", vcat(string.(names), ["distance"]), permutedims(reduce(hcat, [vcat(p, v) for (p, v) in opt.table])))
end

# ---------------------------------------------------------------- validation on held-out conditions
if "validate" in parts
    cal = Dict(String(k) => v for (k, v) in zip(readdlm(joinpath(OUT, "fig7_calibration.csv"), ','; skipstart = 1)[:, 1], readdlm(joinpath(OUT, "fig7_calibration.csv"), ','; skipstart = 1)[:, 2]))
    θ = [cal[string(n)] for n in names]
    println("parameters: ", round.(θ; sigdigits = 3))
    # (a) fate fractions at all three concentrations, several seeds
    rows = Any[]; rows_t = Any[]
    for d in sort(collect(keys(obs)))
        for seed in 1:4
            r = experiment(θ, d; seed = 100 + seed)
            f = fates(r)
            push!(rows, [d, obs[d].label, d in TRAIN ? "train" : "held-out", seed, f.died, f.divided, f.survived, f.n, obs[d].died, obs[d].divided, obs[d].survived, obs[d].n])
            lt = r.lineage
            div = [lt.end_time[i] - lt.birth_time[i] for i in eachindex(lt.id) if lt.fate[i] == :divided && lt.birth_time[i] >= T_FREE]
            dead = [lt.end_time[i] - max(lt.birth_time[i], T_FREE) for i in eachindex(lt.id) if lt.fate[i] == :died && lt.end_time[i] > T_FREE]
            push!(rows_t, [d, seed, "division", length(div), isempty(div) ? NaN : mean(div), isempty(div) ? NaN : std(div)])
            push!(rows_t, [d, seed, "death", length(dead), isempty(dead) ? NaN : mean(dead), isempty(dead) ? NaN : std(dead)])
            @printf("(a) %5.1f µM seed %d: died %.3f divided %.3f survived %.3f (observed %.3f %.3f %.3f)\n", d, seed, f.died, f.divided, f.survived, obs[d].died, obs[d].divided, obs[d].survived)
        end
    end
    save_csv("fig7a_fates.csv", ["cisplatin_uM", "condition", "role", "seed", "died", "divided", "survived", "n_cells", "obs_died", "obs_divided", "obs_survived", "obs_n"], permutedims(reduce(hcat, rows)))
    save_csv("fig7d_timing.csv", ["cisplatin_uM", "seed", "event", "n", "mean_h", "sd_h"], permutedims(reduce(hcat, rows_t)))
    # (b) kill curves over a week at the three concentrations and at the concentration giving 64% death (HCT116-like)
    rows = Any[]
    for d in vcat(sort(collect(keys(obs))), [11.5])
        r = experiment(θ, d; seed = 200, t_drug = 168.0, record_every = 8)
        i0 = argmin(abs.(r.t .- T_FREE))
        for i in eachindex(r.t)
            push!(rows, [d, r.t[i] - T_FREE, r.popsize[i] / r.popsize[i0]])
        end
        f = fates(r)
        @printf("(b) %5.1f µM: N(72h)/N(0) = %.3f, died %.3f, N(168h)/N(0) = %.3f\n", d, r.popsize[argmin(abs.(r.t .- (T_FREE + 72)))] / r.popsize[i0], f.died, r.popsize[end] / r.popsize[i0])
    end
    save_csv("fig7b_killcurves.csv", ["cisplatin_uM", "t_since_drug_h", "N_over_N0"], permutedims(reduce(hcat, rows)))
    # (c) lineage correlations of end fates for sisters, first, second and third cousins (cells present at drug addition)
    function phi(lt, pairs, t_on, t_end)
        a = Float64[]; b = Float64[]
        alive_at(i) = lt.birth_time[i] <= t_on && (isnan(lt.end_time[i]) || lt.end_time[i] > t_on)
        died(i) = lt.fate[i] == :died && lt.end_time[i] <= t_end ? 1.0 : 0.0
        for (i, j) in pairs
            (alive_at(i) && alive_at(j)) || continue
            push!(a, died(i)); push!(b, died(j)); push!(a, died(j)); push!(b, died(i))
        end
        length(a) < 20 && return (NaN, length(a) ÷ 2)
        (std(a) == 0 || std(b) == 0) ? (0.0, length(a) ÷ 2) : (cor(a, b), length(a) ÷ 2)
    end
    rows = Any[]
    for seed in 1:3
        r = experiment(θ, 11.5; N = 60, seed = 300 + seed, t_free = 4.5CYCLE, record_every = 100)   # 4.5 drug-free generations so that third cousins exist
        lt = r.lineage; t_on = 4.5CYCLE; t_end = t_on + T_DRUG
        for g in 1:4
            c, n = phi(lt, kin_pairs(lt, g), t_on, t_end)
            push!(rows, [seed, g, ("sisters", "first cousins", "second cousins", "third cousins")[g], c, n])
        end
        # unrelated pairs (different founders) as the null
        ids = [lt.id[i] for i in eachindex(lt.id) if lt.birth_time[i] <= t_on && (isnan(lt.end_time[i]) || lt.end_time[i] > t_on)]
        rng = Xoshiro(seed); pairs = Tuple{Int,Int}[]
        while length(pairs) < 3000
            i, j = rand(rng, ids), rand(rng, ids)
            lt.clone[i] != lt.clone[j] && push!(pairs, (i, j))
        end
        c, n = phi(lt, pairs, t_on, t_end)
        push!(rows, [seed, 0, "unrelated", c, n])
        @printf("(c) seed %d: %s\n", seed, join((@sprintf("%s %.3f", rows[end-4+k][3], rows[end-4+k][4]) for k in 0:4), ", "))
    end
    save_csv("fig7c_kin_correlation.csv", ["seed", "generations_back", "relation", "fate_correlation", "n_pairs"], permutedims(reduce(hcat, rows)))
end
println("fig7 done")
