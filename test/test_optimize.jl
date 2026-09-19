@testset "schedules, outcomes and optimisation" begin
    # finite cycles
    pd = PulsedDose(1.0; on = 2.0, off = 2.0, start = 0.0, cycles = 2)
    @test dose(pd, 1.0) == 1.0 && dose(pd, 5.0) == 1.0 && dose(pd, 9.0) == 0.0
    # daily boluses and cycle days
    days = cycle_days(5, 28, 2)
    @test days == vcat(0:4, 28:32)
    pk = daily_boluses(2.0, days, log(2) / 2.0; interval = 24.0)
    @test dose(pk, 0.0) ≈ 2.0 && dose(pk, 2.0) ≈ 1.0 && dose(pk, 24.0) ≈ 2.0 + 2.0 * 2.0^(-12)
    @test cumulative_dose(ConstantDose(0.5), 0.0, 10.0) ≈ 5.0 atol = 0.05
    # adaptive therapy: on while above the lower threshold, re-applied above the upper one
    ad = AdaptiveDose(1.0; on_above = 1.0, off_below = 0.5, start = 0.0)
    Ns = [100.0, 80.0, 49.0, 60.0, 101.0, 90.0]
    @test [dose(ad, t, N) for (t, N) in zip(0:5, Ns)] == [1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
    @test dose(ad, 0.0, 100.0) == 1.0           # restarting at an earlier time resets the reference
    @test dose(ConstantDose(2.0), 1.0, 55.0) == 2.0
    # outcomes on a simulated kill-and-regrowth curve
    bd = birth_death_model(k = 0.0, γ = 0.0)
    pert = Perturbation(PulsedDose(1.0; on = 10.0, off = 100.0, cycles = 1); effects = [DeathHazard(h_max = 0.3, EC50 = 1.0, m = 1.0)])
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 10), size_control = AgeTimer(10.0), control = FreeGrowth(max_cells = 5000), threads = false)
    r = simulate_population(bd, [0], 500, (0.0, 40.0); settings = st, perturbation = pert, rng = Xoshiro(1))
    @test log_kill(r) > 0.3
    tp, progressed = time_to_progression(r; threshold = 1.2)
    @test progressed && 10.0 < tp < 40.0
    @test net_growth_rate(r; from = 10.0) > 0.0
    @test cumulative_dose(r) ≈ 10.0 atol = 0.5
    @test extinction_probability([r, r]) == 0.0
    # feedback schedule inside a simulation switches the drug off after the kill
    ad2 = AdaptiveDose(1.0; on_above = 1.0, off_below = 0.5)
    r2 = simulate_population(bd, [0], 500, (0.0, 30.0); settings = st, perturbation = Perturbation(ad2; effects = [DeathHazard(h_max = 0.5, EC50 = 1.0, m = 1.0)]), rng = Xoshiro(2))
    @test any(==(0.0), r2.dose) && r2.dose[1] == 1.0
    # optimiser recovers the minimum of a smooth objective (common random numbers ignore the seed here)
    f(θ, seed) = (log(θ[1]) - log(3.0))^2 + (θ[2] - 0.4)^2
    opt = optimize_schedule(f, [0.1, 0.0], [100.0, 1.0]; names = [:period, :duty], n_grid = 6, log_scale = [true, false], maxiter = 200)
    @test abs(opt.params[1] - 3.0) < 0.3 && abs(opt.params[2] - 0.4) < 0.05
    @test opt.value < 0.01 && length(opt.table) > 36
    @test occursin("period", string(opt))
end
