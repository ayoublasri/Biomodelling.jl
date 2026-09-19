# Figure 8: dose and schedule optimisation against clinical regimens. Time in hours.
#   (a-c) BRAF-mutant melanoma: continuous versus intermittent (SWOG S1320) dosing for pre-existing tolerance
#         without and with a fitness cost of the resistant state, adaptive therapy, and an optimised schedule
#   (d-e) glioblastoma / temozolomide: standard 5/28 versus dose-dense 21/28 regimens with one-compartment
#         pharmacokinetics, for MGMT-methylated and unmethylated populations, and the best regimen at equal cumulative dose
# Usage: julia fig8_schedules.jl [melanoma] [gbm]
include(joinpath(@__DIR__, "common.jl"))
const parts = isempty(ARGS) ? ["melanoma", "gbm"] : ARGS
const HOUR = 1.0; const DAY = 24.0; const WEEK = 7DAY

# ---------------------------------------------------------------- melanoma (WM989-like; Shaffer et al. 2017, S1320)
if "melanoma" in parts
    const CYCLE_M = 4WEEK                                    # net tumour doubling time ≈ 4 weeks (division minus loss)
    const λm = log(2) / CYCLE_M
    const MEMORY_M = 5CYCLE_M                                # memory of the pre-resistant state ≈ 5 generations (MemorySeq)
    const P_ON = 0.005                                       # pre-resistant fraction ≈ 1:200
    k_off_m = (1 - P_ON) / MEMORY_M; k_on_m = P_ON / MEMORY_M
    mel = telegraph_model(k_on = k_on_m, k_off = k_off_m, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2, gene = :R, mrna = :mRNA, protein = :P)
    xm = initial_state(mel; R_off = 1)
    # drug effects at dose 1 (fractional-killing dose): sensitive cells arrest and die with a half-life of about two weeks,
    # resistant (P-high) cells keep growing; the drug stabilises the resistant state ten-fold (drug-induced reprogramming)
    death_m = DeathHazard(h_max = 0.0023, EC50 = 0.3, m = 2.0, protect = :P, K = 150.0, q = 4.0)
    arrest_m = GrowthInhibition(IC50 = 0.3, m = 2.0, protect = :P, K = 150.0, q = 4.0)
    stabilise = RateModulation(:k_off, d -> 1 / (1 + 9d))
    cost = GrowthCost(:P; K = 150.0, q = 4.0, max_cost = 0.5)  # drug addiction / fitness cost of resistant cells (Das Thakur et al. 2013)
    mechanisms = [("no fitness cost", DrugEffect[death_m, arrest_m, stabilise]), ("fitness cost", DrugEffect[death_m, arrest_m, stabilise, cost])]
    burn_m = PopulationSettings(dt = 4.0, growth = ExponentialGrowth(λm), size_control = Sizer(2.0; cv = 0.1), control = ConstantN(), track_lineage = false, record_every = 10_000)
    st_m = PopulationSettings(dt = 4.0, growth = ExponentialGrowth(λm), size_control = Sizer(2.0; cv = 0.1), control = FreeGrowth(max_cells = 20_000), track_lineage = false, record_every = 6)
    const N_M = 2000; const LEAD_IN = 8WEEK; const T_M = 32WEEK
    founders = let b = simulate_population(mel, xm, N_M, (0.0, 10CYCLE_M); settings = burn_m, rng = Xoshiro(11)); final_snapshot(b) end
    function treat_m(sched, effects; seed = 1, T = T_M)
        simulate_population(mel, founders.counts, N_M, (0.0, T); settings = st_m, V0 = founders.volume, perturbation = Perturbation(sched; effects), rng = Xoshiro(seed))
    end
    ttp_m(r) = time_to_progression(r; from = LEAD_IN, threshold = 1.73)      # 20 % diameter increase from the nadir ≈ 1.73-fold in cell number, after randomisation
    # clinical schedules: continuous; S1320 intermittent = 8-week continuous lead-in, then 3 weeks off / 5 weeks on
    continuous = ConstantDose(1.0)
    s1320 = PiecewiseDose([0.0, 8WEEK, 11WEEK, 16WEEK, 19WEEK, 24WEEK, 27WEEK], [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
    schedules = [("continuous", () -> continuous), ("intermittent (S1320)", () -> s1320), ("adaptive (50 %)", () -> AdaptiveDose(1.0; on_above = 1.0, off_below = 0.5))]
    rows = Any[]; rows_t = Any[]
    for (mtag, eff) in mechanisms, (stag, mk) in schedules, seed in 1:4
        r = treat_m(mk(), eff; seed = 20 + seed)
        ttp, prog = ttp_m(r)
        push!(rows, [mtag, stag, seed, (ttp - LEAD_IN) / WEEK, prog, r.popsize[end] / r.popsize[1], log_kill(r), cumulative_dose(r) / WEEK, minimum(r.popsize) / r.popsize[1]])
        seed == 1 && for i in eachindex(r.t); push!(rows_t, [mtag, stag, r.t[i] / WEEK, r.popsize[i] / r.popsize[1], r.dose[i]]); end
        @printf("(a) %-16s %-22s seed %d  TTP after randomisation %.1f weeks%s  N_end/N0 %.3f  nadir %.3f  cumulative dose %.1f weeks\n", mtag, stag, seed, (ttp - LEAD_IN) / WEEK, prog ? "" : " (censored)", r.popsize[end] / r.popsize[1], minimum(r.popsize) / r.popsize[1], cumulative_dose(r) / WEEK)
    end
    save_csv("fig8a_melanoma_schedules.csv", ["mechanism", "schedule", "seed", "ttp_weeks", "progressed", "N_end_over_N0", "log_kill", "cumulative_dose_weeks", "nadir_over_N0"], permutedims(reduce(hcat, rows)))
    save_csv("fig8b_melanoma_trajectories.csv", ["mechanism", "schedule", "t_weeks", "N_over_N0", "dose"], permutedims(reduce(hcat, rows_t)))
    # (c) optimised intermittent schedule (period, duty cycle) after the lead-in, for each mechanism
    for (mtag, eff) in mechanisms
        function objective(θ, seed)
            period, duty = θ
            pulses = PulsedDose(1.0; on = duty * period, off = (1 - duty) * period, start = LEAD_IN)
            r = treat_m(FunctionDose(t -> t < LEAD_IN ? 1.0 : dose(pulses, t)), eff; seed = 40 + seed)
            -(ttp_m(r)[1] - LEAD_IN) / WEEK
        end
        opt = optimize_schedule(objective, [1WEEK, 0.2], [12WEEK, 1.0]; names = [:period_h, :duty], n_grid = 4, seeds = 1:1, maxiter = 16, log_scale = [true, false], rng = Xoshiro(5))
        @printf("(c) %-16s best period %.1f weeks, duty %.2f → TTP %.1f weeks\n", mtag, opt.params[1] / WEEK, opt.params[2], -opt.value)
        save_csv("fig8c_melanoma_optimum_$(replace(mtag, " " => "_")).csv", ["period_weeks", "duty", "ttp_weeks"], permutedims(reduce(hcat, [[p[1] / WEEK, p[2], -v] for (p, v) in opt.table])))
    end
end

# ---------------------------------------------------------------- glioblastoma / temozolomide (N15-0385-like; RTOG 0525)
if "gbm" in parts
    const CYCLE_G = 40DAY                                    # net tumour doubling time ≈ 40 days
    const λg = log(2) / CYCLE_G
    const K_E = log(2) / 2.1                                 # temozolomide elimination half-life 2.1 h (Ostermann et al. 2004)
    memory_g = 4CYCLE_G
    mgmt(p_on) = telegraph_model(k_on = p_on / memory_g, k_off = (1 - p_on) / memory_g, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2, gene = :MGMT, mrna = :mRNA, protein = :P)
    death_g = DeathHazard(h_max = 0.03, EC50 = 0.4, m = 2.0, protect = :P, K = 150.0, q = 4.0)   # dose in units of the standard daily bolus peak
    burn_g = PopulationSettings(dt = 4.0, growth = ExponentialGrowth(λg), size_control = Sizer(2.0; cv = 0.1), control = ConstantN(), track_lineage = false, record_every = 10_000)
    st_g = PopulationSettings(dt = 1.0, growth = ExponentialGrowth(λg), size_control = Sizer(2.0; cv = 0.1), control = FreeGrowth(max_cells = 20_000), track_lineage = false, record_every = 24)
    const N_G = 1500; const N_CYCLES = 6; const T_G = N_CYCLES * 28DAY
    regimens = [("standard 5/28", daily_boluses(1.0, cycle_days(5, 28, N_CYCLES), K_E)),                      # 200 mg/m² days 1-5
                ("dose-dense 21/28", daily_boluses(0.5, cycle_days(21, 28, N_CYCLES), K_E)),                  # 100 mg/m² days 1-21 (2.1× cumulative)
                ("dense, equal cumulative", daily_boluses(5 / 21, cycle_days(21, 28, N_CYCLES), K_E))]         # 21 days at the standard cumulative dose
    rows = Any[]; rows_t = Any[]
    for (ptag, p_on) in (("MGMT methylated (1 % expressing)", 0.01), ("MGMT unmethylated (30 % expressing)", 0.30))
        m = mgmt(p_on); x0 = initial_state(m; MGMT_off = 1)
        founders_g = final_snapshot(simulate_population(m, x0, N_G, (0.0, 8CYCLE_G); settings = burn_g, rng = Xoshiro(13)))
        run_g(sched; seed) = simulate_population(m, founders_g.counts, N_G, (0.0, T_G); settings = st_g, V0 = founders_g.volume, perturbation = Perturbation(sched; effects = [death_g]), rng = Xoshiro(seed))
        for (rtag, sched) in regimens, seed in 1:2
            r = run_g(sched; seed = 60 + seed)
            push!(rows, [ptag, rtag, seed, r.popsize[end] / r.popsize[1], log_kill(r), net_growth_rate(r) * WEEK, cumulative_dose(sched, 0.0, T_G) / DAY, time_to_progression(r; threshold = 1.73)[1] / WEEK])
            seed == 1 && for i in eachindex(r.t); push!(rows_t, [ptag, rtag, r.t[i] / WEEK, r.popsize[i] / r.popsize[1], r.dose[i]]); end
            @printf("(d) %-36s %-24s seed %d  N_end/N0 %.3f  log kill %.2f  cumulative %.1f bolus-days\n", ptag, rtag, seed, r.popsize[end] / r.popsize[1], log_kill(r), cumulative_dose(sched, 0.0, T_G) / DAY)
        end
        # (e) number of dosing days per 28-day cycle at the standard cumulative dose (5 bolus units per cycle)
        function objective(θ, seed)
            days_on = clamp(round(Int, θ[1]), 1, 28)
            sched = daily_boluses(5 / days_on, cycle_days(days_on, 28, N_CYCLES), K_E)
            r = run_g(sched; seed = 80 + seed)
            log(max(r.popsize[end], 1.0) / r.popsize[1])
        end
        opt = optimize_schedule(objective, [1.0], [28.0]; names = [:days_on], n_grid = 8, seeds = 1:1, refine = false, log_scale = false)
        @printf("(e) %-36s best days on per cycle ≈ %.0f (log N_end/N0 %.2f)\n", ptag, opt.params[1], opt.value)
        save_csv("fig8e_gbm_days_on_$(startswith(ptag, "MGMT m") ? "methylated" : "unmethylated").csv", ["days_on", "log_N_end_over_N0"], permutedims(reduce(hcat, [[round(p[1]), v] for (p, v) in opt.table])))
    end
    save_csv("fig8d_gbm_regimens.csv", ["population", "regimen", "seed", "N_end_over_N0", "log_kill", "net_growth_per_week", "cumulative_bolus_days", "ttp_weeks"], permutedims(reduce(hcat, rows)))
    save_csv("fig8d_gbm_trajectories.csv", ["population", "regimen", "t_weeks", "N_over_N0", "dose"], permutedims(reduce(hcat, rows_t)))
end
println("fig8 done")
