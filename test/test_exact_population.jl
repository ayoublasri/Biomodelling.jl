@testset "exact population and lineage laws" begin
    # A model whose stationary distribution is known in closed form in both the single-lineage
    # and the population-snapshot setting: zero-order production, first-order decay, and a
    # memoryless interdivision time with binomial partitioning at division. Writing the
    # generating function as G(u) = sum_j c_j u^j with u = z - 1, the factorial moments obey
    #   lineage     c_j = k c_{j-1} / (gamma j + lambda (1 - 2^-j))
    #   population  c_j = k c_{j-1} / (gamma j + lambda (2 - 2^(1-j)))
    # The two differ because a snapshot of a growing population over-weights cells that have
    # just divided, and so have just lost half their molecules.
    k, γ, λd = 20.0, 1.0, 1.0
    facmom(den) = (c = [1.0]; for j in 1:2; push!(c, k * c[end] / den(j)); end; c)
    lin = facmom(j -> γ * j + λd * (1 - 2.0^-j))
    pop = facmom(j -> γ * j + λd * (2 - 2.0^(1 - j)))
    mean_var(c) = (c[2], 2c[3] + c[2] - c[2]^2)
    m_lin, v_lin = mean_var(lin)
    m_pop, v_pop = mean_var(pop)
    @test m_lin ≈ 40 / 3
    @test m_pop ≈ 10.0

    m = ReactionModel([Reaction("production", [], [:M], MassAction(:k); volume = :none),
                       Reaction("decay", [:M], [], MassAction(:g); volume = :none)];
                      params = (k = k, g = γ))
    settings(control) = PopulationSettings(; dt = 0.005, kernel = DirectSSA(), growth = ExponentialGrowth(0.0),
                                           size_control = AgeTimer(1 / λd; cv = 1.0, dist = :gamma),
                                           partitioning = BinomialPartition(), control = control,
                                           track_lineage = false, record_every = 10_000, threads = false)

    # mother machine: one daughter kept at random at every division, so the cells are
    # independent lineages and the sample follows the single-lineage law
    rl = simulate_population(m, initial_state(m), 3000, (0.0, 8.0); settings = settings(MotherMachine()), rng = Xoshiro(1))
    @test all(size(c, 1) == 3000 for c in rl.counts)
    xl = final_snapshot(rl).counts[:, 1]
    @test abs(mean(xl) - m_lin) < 0.5
    @test abs(std(xl) - sqrt(v_lin)) < 0.5

    # freely growing population: the same cells, sampled as a snapshot
    rp = simulate_population(m, initial_state(m), 200, (0.0, 8.0); settings = settings(FreeGrowth(max_cells = 4000)), rng = Xoshiro(2))
    xp = final_snapshot(rp).counts[:, 1]
    @test abs(mean(xp) - m_pop) < 1.0
    @test mean(xp) < mean(xl) - 1.5                      # the two modes are distinguishable

    # gamma-distributed interdivision times have the requested mean and spread
    for cv in (0.25, 1.0)
        ts = [BM.sample_target(AgeTimer(3.0; cv = cv, dist = :gamma), 1.0, Xoshiro(i)) for i in 1:20_000]
        @test abs(mean(ts) - 3.0) < 0.15
        @test abs(std(ts) / mean(ts) - cv) < 0.05
    end

    # founder ages are taken from age0, and a deterministic timer keeps them on that grid
    tm = telegraph_model(k_on = 1.0, k_off = 1.0, k_tx = 20.0, k_dm = 1.0, k_tl = 2.0, k_dp = 0.1)
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 10), size_control = AgeTimer(10.0),
                            control = MotherMachine(), track_lineage = false, threads = false)
    ages = [0.0, 2.5, 5.0, 7.5]
    r0 = simulate_population(tm, initial_state(tm; G_off = 1), 4, (0.0, 0.0); settings = st,
                             rng = Xoshiro(3), age0 = ages)
    @test r0.age[1] ≈ ages
    r1 = simulate_population(tm, initial_state(tm; G_off = 1), 4, (0.0, 30.0); settings = st,
                             rng = Xoshiro(3), age0 = ages)
    @test sort(unique(round.(r1.age[end]; digits = 6))) ⊆ [0.0, 2.5, 5.0, 7.5]
end
