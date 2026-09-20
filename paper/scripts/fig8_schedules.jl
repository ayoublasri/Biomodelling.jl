# Figure 8: dose and schedule optimisation against clinical regimens. Time in hours.
#   (a-c) BRAF-mutant melanoma: continuous versus intermittent (SWOG S1320) dosing for pre-existing tolerance
#         without and with a fitness cost of the resistant state, adaptive therapy, and an optimised schedule
#   (d-e) glioblastoma / temozolomide: standard 5/28 versus dose-dense 21/28 regimens with one-compartment
#         pharmacokinetics, for MGMT-methylated and unmethylated populations, and the best regimen at equal cumulative dose
# Usage: julia fig8_schedules.jl [melanoma] [gbm]
include(joinpath(@__DIR__, "common.jl"))
const parts = isempty(ARGS) ? ["melanoma", "gbm"] : ARGS
const HOUR = 1.0; const DAY = 24.0; const WEEK = 7DAY
"""Founders with promoter states drawn from the stationary distribution (memory is long, so a burn-in from a
common state would not reach it); a short burn-in then randomises volumes and ages."""
function stationary_founders(m, N, p_on, burn, seed; cycle)
    rng = Xoshiro(seed); X = zeros(Int, N, nspecies(m))
    ioff, ion = m.promoter_groups[1][1], m.promoter_groups[1][2]        # promoter groups hold species indices, [off, on]
    imr = speciesindex(m, :mRNA); ip = speciesindex(m, :P)
    for i in 1:N
        if rand(rng) < p_on
            X[i, ion] = 1; X[i, imr] = 30; X[i, ip] = 600
        else
            X[i, ioff] = 1
        end
    end
    final_snapshot(simulate_population(m, X, N, (0.0, 3cycle); settings = burn, rng = Xoshiro(seed + 1)))
end


# ---------------------------------------------------------------- melanoma (WM989-like; Shaffer et al. 2017, S1320)
if ("melanoma" in parts) || ("melanoma_cycle" in parts)
    const CYCLE_M = 4WEEK                                    # net tumour doubling time ≈ 4 weeks (division minus loss)
    const λm = log(2) / CYCLE_M
    const MEMORY_M = 5CYCLE_M                                # memory of the pre-resistant state ≈ 5 generations (MemorySeq)
    const P_ON = 0.005                                       # pre-resistant fraction ≈ 1:200
    k_off_m = (1 - P_ON) / MEMORY_M; k_on_m = P_ON / MEMORY_M
    mel = telegraph_model(k_on = k_on_m, k_off = k_off_m, k_tx = 3.0, k_dm = 0.1, k_tl = 0.4, k_dp = 0.02, gene = :R, mrna = :mRNA, protein = :P)
    xm = initial_state(mel; R_off = 1)
    # drug effects at dose 1 (fractional-killing dose): sensitive cells arrest and die with a half-life of about two weeks,
    # resistant (P-high) cells keep growing; the drug stabilises the resistant state ten-fold (drug-induced reprogramming)
    death_m = DeathHazard(h_max = 0.0023, EC50 = 0.3, m = 2.0, protect = :P, K = 150.0, q = 4.0)
    arrest_m = GrowthInhibition(IC50 = 0.3, m = 2.0, protect = :P, K = 150.0, q = 4.0)
    stabilise = RateModulation(:k_off, d -> 1 / (1 + 9d))
    cost = GrowthCost(:P; K = 150.0, q = 4.0, max_cost = 0.5)  # drug addiction / fitness cost of resistant cells (Das Thakur et al. 2013)
    slow_all = GrowthInhibition(IC50 = 1.0, m = 2.0)             # partial protection: resistant cells still grow at half speed under drug
    mechanisms = [("no fitness cost", DrugEffect[death_m, arrest_m, stabilise]), ("fitness cost", DrugEffect[death_m, arrest_m, stabilise, cost]),
                  ("partial protection", DrugEffect[death_m, arrest_m, slow_all, stabilise])]
    burn_m = PopulationSettings(dt = 4.0, growth = ExponentialGrowth(λm), size_control = Sizer(2.0; cv = 0.1), control = ConstantN(), track_lineage = false, record_every = 10_000)
    st_m = PopulationSettings(dt = 4.0, growth = ExponentialGrowth(λm), size_control = Sizer(2.0; cv = 0.1), control = FreeGrowth(max_cells = 20_000), track_lineage = false, record_every = 6)
    const N_M = 2000; const LEAD_IN = 8WEEK; const T_M = 60WEEK
    founders = stationary_founders(mel, N_M, P_ON, burn_m, 11; cycle = CYCLE_M)
    function treat_m(sched, effects; seed = 1, T = T_M)
        simulate_population(mel, founders.counts, N_M, (0.0, T); settings = st_m, V0 = founders.volume, perturbation = Perturbation(sched; effects), rng = Xoshiro(seed))
    end
    ttp_m(r) = time_to_progression(r; from = LEAD_IN, threshold = 1.73)      # 20 % diameter increase from the nadir ≈ 1.73-fold in cell number, after randomisation
    """Time after randomisation at which the population first exceeds 120 % of its pre-treatment size (loss of control),
    the end point used for adaptive therapy; censored at the end of the simulation."""
    function ttp_baseline(r)
        i = findfirst(k -> r.t[k] >= LEAD_IN && r.popsize[k] >= 1.2 * r.popsize[1], eachindex(r.t))
        i === nothing ? (r.t[end], false) : (r.t[i], true)
    end
    # clinical schedules: continuous; S1320 intermittent = 8-week continuous lead-in, then 3 weeks off / 5 weeks on
    continuous = ConstantDose(1.0)
    s1320 = let ts = [0.0, LEAD_IN]; ds = [1.0, 0.0]; t = LEAD_IN
        while t < T_M                                            # 3 weeks off, 5 weeks on, repeated
            t += 3WEEK; push!(ts, t); push!(ds, 1.0); t += 5WEEK; push!(ts, t); push!(ds, 0.0)
        end
        PiecewiseDose(ts, ds)
    end
    schedules = [("continuous", () -> continuous), ("intermittent (S1320)", () -> s1320), ("adaptive (50 %)", () -> AdaptiveDose(1.0; on_above = 1.0, off_below = 0.5))]
if "melanoma" in parts
    rows = Any[]; rows_t = Any[]
    for (mtag, eff) in mechanisms, (stag, mk) in schedules, seed in 1:4
        r = treat_m(mk(), eff; seed = 20 + seed)
        ttp, prog = ttp_m(r); tb, progb = ttp_baseline(r)
        push!(rows, [mtag, stag, seed, (ttp - LEAD_IN) / WEEK, prog, (tb - LEAD_IN) / WEEK, progb, r.popsize[end] / r.popsize[1], log_kill(r), cumulative_dose(r) / WEEK, minimum(r.popsize) / r.popsize[1]])
        seed == 1 && for i in eachindex(r.t); push!(rows_t, [mtag, stag, r.t[i] / WEEK, r.popsize[i] / r.popsize[1], r.dose[i]]); end
        @printf("(a) %-18s %-22s seed %d  TTP nadir %.1f w%s  loss of control %.1f w%s  N_end/N0 %.3f  nadir %.3f  dose %.1f w\n", mtag, stag, seed, (ttp - LEAD_IN) / WEEK, prog ? "" : "*", (tb - LEAD_IN) / WEEK, progb ? "" : "*", r.popsize[end] / r.popsize[1], minimum(r.popsize) / r.popsize[1], cumulative_dose(r) / WEEK)
    end
    save_csv("fig8a_melanoma_schedules.csv", ["mechanism", "schedule", "seed", "ttp_nadir_weeks", "progressed_nadir", "ttp_baseline_weeks", "progressed_baseline", "N_end_over_N0", "log_kill", "cumulative_dose_weeks", "nadir_over_N0"], permutedims(reduce(hcat, rows)))
    save_csv("fig8b_melanoma_trajectories.csv", ["mechanism", "schedule", "t_weeks", "N_over_N0", "dose"], permutedims(reduce(hcat, rows_t)))
    # (c) optimised intermittent schedule (period, duty cycle) after the lead-in, for each mechanism
    for (mtag, eff) in mechanisms
        function objective(θ, seed)
            period, duty = θ
            pulses = PulsedDose(1.0; on = duty * period, off = (1 - duty) * period, start = LEAD_IN)
            r = treat_m(FunctionDose(t -> t < LEAD_IN ? 1.0 : dose(pulses, t)), eff; seed = 40 + seed)
            -(ttp_baseline(r)[1] - LEAD_IN) / WEEK
        end
        opt = optimize_schedule(objective, [1WEEK, 0.2], [12WEEK, 1.0]; names = [:period_h, :duty], n_grid = 4, seeds = 1:1, maxiter = 16, log_scale = [true, false], rng = Xoshiro(5))
        @printf("(c) %-18s best period %.1f weeks, duty %.2f → loss of control after %.1f weeks\n", mtag, opt.params[1] / WEEK, opt.params[2], -opt.value)
        save_csv("fig8c_melanoma_optimum_$(replace(mtag, " " => "_")).csv", ["period_weeks", "duty", "ttp_weeks"], permutedims(reduce(hcat, [[p[1] / WEEK, p[2], -v] for (p, v) in opt.table])))
    end
end  # "melanoma"

    # Robustness of the schedule ranking to cell-cycle-dependent killing. A BRAF/MEK inhibitor does not make
    # replication-coupled lesions as a platinum drug does, but it does act on cells traversing the cycle, so this
    # asks whether the ranking survives a hazard of which only a fraction acts outside a window of the cycle.
    if "melanoma_cycle" in parts
        CYC_M = CycleSensitivity(baseline = 0.25, center = 0.5, width = 0.15)
        rows_y = Any[]
        for (mtag, eff) in mechanisms, (stag, mk) in schedules, seed in 1:4
            r = treat_m(mk(), vcat(eff, DrugEffect[CYC_M]); seed = 20 + seed)
            tb, progb = ttp_baseline(r)
            push!(rows_y, [mtag, stag, seed, (tb - LEAD_IN) / WEEK, progb, r.popsize[end] / r.popsize[1], cumulative_dose(r) / WEEK])
            @printf("(cycle) %-18s %-22s seed %d  loss of control %.1f w%s\n", mtag, stag, seed, (tb - LEAD_IN) / WEEK, progb ? "" : "*")
        end
        save_csv("fig8f_melanoma_cycle.csv", ["mechanism", "schedule", "seed", "ttp_baseline_weeks", "progressed_baseline", "N_end_over_N0", "cumulative_dose_weeks"], permutedims(reduce(hcat, rows_y)))
    end
end

# ---------------------------------------------------------------- glioblastoma / temozolomide (N15-0385-like; RTOG 0525)
if "gbm" in parts
    const CYCLE_G = 40DAY                                    # net tumour doubling time ≈ 40 days
    const λg = log(2) / CYCLE_G
    const K_E = log(2) / 2.1                                 # temozolomide elimination half-life 2.1 h (Ostermann et al. 2004)
    memory_g = 4CYCLE_G
    mgmt(p_on) = telegraph_model(k_on = p_on / memory_g, k_off = (1 - p_on) / memory_g, k_tx = 3.0, k_dm = 0.1, k_tl = 0.4, k_dp = 0.02, gene = :MGMT, mrna = :mRNA, protein = :P)
    death_g = DeathHazard(h_max = 0.03, EC50 = 0.4, m = 2.0, protect = :P, K = 150.0, q = 4.0)   # dose in units of the standard daily bolus peak
    burn_g = PopulationSettings(dt = 4.0, growth = ExponentialGrowth(λg), size_control = Sizer(2.0; cv = 0.1), control = ConstantN(), track_lineage = false, record_every = 10_000)
    st_g = PopulationSettings(dt = 1.0, growth = ExponentialGrowth(λg), size_control = Sizer(2.0; cv = 0.1), control = FreeGrowth(max_cells = 20_000), track_lineage = false, record_every = 24)
    const N_G = 1500; const N_CYCLES = 6; const T_G = N_CYCLES * 28DAY
    regimens = [("standard 5/28", daily_boluses(1.0, cycle_days(5, 28, N_CYCLES), K_E)),                      # 200 mg/m² days 1-5
                ("dose-dense 21/28", daily_boluses(0.5, cycle_days(21, 28, N_CYCLES), K_E)),                  # 100 mg/m² days 1-21 (2.1× cumulative)
                ("dense, equal cumulative", daily_boluses(5 / 21, cycle_days(21, 28, N_CYCLES), K_E))]         # 21 days at the standard cumulative dose
    # MGMT is a suicide enzyme: one molecule is spent per lesion repaired, and lesions form in proportion to
    # the dose, so the pool follows the cumulative exposure rather than the peak concentration. k is set so the
    # consumption rate at the standard bolus peak matches the first-order parameterisation it replaces
    # (k_dp -> k_dp(1 + 20d)), which isolates the difference in how consumption integrates over time.
    deplete = SuicideConsumption(:P; k = 300.0, K_m = 150.0)
    mechanisms_g = [("MGMT stable", DrugEffect[death_g]), ("MGMT consumed by drug", DrugEffect[death_g, deplete])]
    rows = Any[]; rows_t = Any[]
    for (ptag, p_on) in (("MGMT methylated (1 % expressing)", 0.01), ("MGMT unmethylated (30 % expressing)", 0.30)), (gtag, eff) in mechanisms_g
        m = mgmt(p_on)
        founders_g = stationary_founders(m, N_G, p_on, burn_g, 13; cycle = CYCLE_G)
        run_g(sched; seed) = simulate_population(m, founders_g.counts, N_G, (0.0, T_G); settings = st_g, V0 = founders_g.volume, perturbation = Perturbation(sched; effects = eff), rng = Xoshiro(seed))
        for (rtag, sched) in regimens, seed in 1:2
            r = run_g(sched; seed = 60 + seed)
            push!(rows, [ptag, gtag, rtag, seed, r.popsize[end] / r.popsize[1], log_kill(r), net_growth_rate(r) * WEEK, cumulative_dose(sched, 0.0, T_G) / DAY, time_to_progression(r; threshold = 1.73)[1] / WEEK])
            seed == 1 && for i in eachindex(r.t); push!(rows_t, [ptag, gtag, rtag, r.t[i] / WEEK, r.popsize[i] / r.popsize[1], r.dose[i]]); end
            @printf("(d) %-36s %-22s %-24s seed %d  N_end/N0 %.3f  log kill %.2f  cumulative %.1f bolus-days\n", ptag, gtag, rtag, seed, r.popsize[end] / r.popsize[1], log_kill(r), cumulative_dose(sched, 0.0, T_G) / DAY)
        end
        # (e) number of dosing days per 28-day cycle at the standard cumulative dose (5 bolus units per cycle)
        function objective(θ, seed)
            days_on = clamp(round(Int, θ[1]), 1, 28)
            sched = daily_boluses(5 / days_on, cycle_days(days_on, 28, N_CYCLES), K_E)
            r = run_g(sched; seed = 80 + seed)
            log(max(r.popsize[end], 1.0) / r.popsize[1])
        end
        opt = optimize_schedule(objective, [1.0], [28.0]; names = [:days_on], n_grid = 8, seeds = 1:1, refine = false, log_scale = false)
        @printf("(e) %-36s %-22s best days on per cycle ≈ %.0f (log N_end/N0 %.2f)\n", ptag, gtag, opt.params[1], opt.value)
        save_csv("fig8e_gbm_days_on_$(startswith(ptag, "MGMT m") ? "methylated" : "unmethylated")_$(startswith(gtag, "MGMT s") ? "stable" : "consumed").csv", ["days_on", "log_N_end_over_N0"], permutedims(reduce(hcat, [[round(p[1]), v] for (p, v) in opt.table])))
    end
    save_csv("fig8d_gbm_regimens.csv", ["population", "mechanism", "regimen", "seed", "N_end_over_N0", "log_kill", "net_growth_per_week", "cumulative_bolus_days", "ttp_weeks"], permutedims(reduce(hcat, rows)))
    save_csv("fig8d_gbm_trajectories.csv", ["population", "mechanism", "regimen", "t_weeks", "N_over_N0", "dose"], permutedims(reduce(hcat, rows_t)))
end
println("fig8 done")
