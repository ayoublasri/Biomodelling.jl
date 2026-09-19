@testset "perturbation layer" begin
    @test dose(ConstantDose(2.0), 5.0) == 2.0
    pd = PulsedDose(1.0; on = 2.0, off = 3.0, start = 1.0)
    @test dose(pd, 0.5) == 0.0 && dose(pd, 1.5) == 1.0 && dose(pd, 3.5) == 0.0 && dose(pd, 6.5) == 1.0
    pw = PiecewiseDose([0.0, 10.0, 20.0], [0.0, 1.0, 0.5])
    @test dose(pw, -1.0) == 0.0 && dose(pw, 5.0) == 0.0 && dose(pw, 15.0) == 1.0 && dose(pw, 25.0) == 0.5
    pk = BolusPK([0.0, 10.0], [1.0, 1.0], 0.1)
    @test dose(pk, 10.0) ≈ exp(-1.0) + 1.0
    @test_throws ArgumentError PiecewiseDose([1.0, 0.0], [1.0, 1.0])
    # death hazard: survival of non-growing cells ≈ exp(-h t)
    bd = birth_death_model(k = 0.0, γ = 0.0)
    pert = Perturbation(ConstantDose(1.0); effects = [DeathHazard(h_max = 0.2, EC50 = 1.0, m = 1.0)])
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(0.0), size_control = AgeTimer(1e9), control = FreeGrowth(), threads = false)
    res = simulate_population(bd, [0], 2000, (0.0, 10.0); settings = st, perturbation = pert, rng = Xoshiro(1))
    @test abs(res.popsize[end] / 2000 - exp(-0.1 * 10)) < 0.05
    @test count(==(:died), res.lineage.fate) == 2000 - round(Int, res.popsize[end])
    # protection by a resistance protein: high-expressing cells survive preferentially
    tm = telegraph_model(k_on = 0.02, k_off = 0.02, k_tx = 20.0, k_dm = 1.0, k_tl = 5.0, k_dp = 0.2)
    x0 = initial_state(tm; G_off = 1)
    pert2 = Perturbation(PiecewiseDose([0.0, 40.0], [0.0, 1.0]);
                         effects = [DeathHazard(h_max = 0.5, EC50 = 0.5, protect = :protein, K = 100.0, q = 4.0)])
    st2 = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20), control = FreeGrowth(max_cells = 4000), threads = false)
    res2 = simulate_population(tm, x0, 300, (0.0, 70.0); settings = st2, perturbation = pert2, rng = Xoshiro(2))
    before = snapshot(res2; t = 39.0); after = final_snapshot(res2)
    pi = speciesindex(tm, :protein)
    @test mean(after.counts[:, pi] ./ after.volume) > 1.5 * mean(before.counts[:, pi] ./ before.volume)
    @test count(==(:died), res2.lineage.fate) > 50
    # growth inhibition slows the population
    pert3 = Perturbation(ConstantDose(1.0); effects = [GrowthInhibition(IC50 = 1.0, m = 1.0)])
    r_ctrl = simulate_population(tm, x0, 100, (0.0, 30.0); settings = st2, rng = Xoshiro(3))
    r_inh = simulate_population(tm, x0, 100, (0.0, 30.0); settings = st2, perturbation = pert3, rng = Xoshiro(3))
    @test r_inh.popsize[end] < 0.8 * r_ctrl.popsize[end]
    # rate modulation and gene knockdown change expression
    pert4 = Perturbation(ConstantDose(1.0); effects = [RateModulation(:k_tx, d -> 1 + 2d)])
    r_mod = simulate_population(tm, x0, 100, (0.0, 40.0); settings = st2, perturbation = pert4, rng = Xoshiro(4))
    r_ref = simulate_population(tm, x0, 100, (0.0, 40.0); settings = st2, rng = Xoshiro(4))
    mi = speciesindex(tm, :mRNA)
    @test mean(r_mod.counts[end][:, mi]) > 2 * mean(r_ref.counts[end][:, mi])
    pert5 = Perturbation(; gene_perturbations = [GenePerturbation(:k_tx, 0.0; fraction = 0.5, t_start = 0.0)])
    r_ko = simulate_population(tm, x0, 400, (0.0, 30.0); settings = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20), threads = false), perturbation = pert5, rng = Xoshiro(5))
    sn = final_snapshot(r_ko)
    @test 0.3 < mean(sn.perturbed) < 0.7
    @test mean(sn.counts[sn.perturbed, mi]) < 0.3 * mean(sn.counts[.!sn.perturbed, mi])
end

@testset "growth cost and per-cell initial states" begin
    tm = telegraph_model(k_on = 0.02, k_off = 0.02, k_tx = 20.0, k_dm = 1.0, k_tl = 5.0, k_dp = 0.2)
    x0 = initial_state(tm; G_off = 1)
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20), control = FreeGrowth(max_cells = 4000), threads = false)
    burn = simulate_population(tm, x0, 200, (0.0, 40.0); settings = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20), threads = false), rng = Xoshiro(1))
    sn = final_snapshot(burn)
    r0 = simulate_population(tm, sn.counts, 200, (0.0, 30.0); settings = st, V0 = sn.volume, rng = Xoshiro(2))
    @test r0.counts[1] == sn.counts && r0.volume[1] == sn.volume
    @test_throws ArgumentError simulate_population(tm, sn.counts[1:10, :], 200, (0.0, 1.0); settings = st)
    cost = Perturbation(; effects = [GrowthCost(:protein; K = 10.0, q = 4.0, max_cost = 0.9)])
    r1 = simulate_population(tm, sn.counts, 200, (0.0, 30.0); settings = st, V0 = sn.volume, perturbation = cost, rng = Xoshiro(2))
    @test r1.popsize[end] < 0.8 * r0.popsize[end]
end
