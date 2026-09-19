@testset "population: growth, division, partitioning" begin
    tm = telegraph_model(k_on = 1.0, k_off = 1.0, k_tx = 40.0, k_dm = 1.0)
    x0 = initial_state(tm; G_off = 1)
    λ = log(2) / 10
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05),
                            partitioning = BinomialPartition(σ = 0.02), replication = Replication(0.5), threads = false)
    res = simulate_population(tm, x0, 200, (0.0, 100.0); settings = st, rng = Xoshiro(1))
    @test length(res.t) == 1001
    @test all(size(c, 1) == 200 for c in res.counts)          # ConstantN
    @test all(all(v .> 0) for v in res.volume)
    lt = res.lineage
    @test all(lt.id .== 1:length(lt))                         # row == id invariant
    @test all(lt.birth_time[lt.parent .!= 0] .>= 0)
    @test all(i -> lt.parent[i] == 0 || lt.birth_time[i] >= lt.birth_time[lt.parent[i]], eachindex(lt.id))
    # cell-cycle time ≈ ln2/λ for a sizer with 2× growth
    cyc = [lt.end_time[i] - lt.birth_time[i] for i in eachindex(lt.id) if lt.fate[i] == :divided && lt.generation[i] > 0]
    @test abs(mean(cyc) - log(2) / λ) < 1.0
    # volumes at division ≈ 2, at birth ≈ 1
    @test abs(mean(lt.V_end[lt.fate .== :divided]) - 2.0) < 0.1
    @test abs(mean(lt.V_birth[lt.parent .!= 0]) - 1.0) < 0.1
    # promoter inheritance: every cell always carries exactly `copies` promoter copies
    for i in (1, 500, 1001)
        sn = snapshot(res, i)
        @test all(sn.counts[:, 1] .+ sn.counts[:, 2] .== sn.copies)
    end
    @test Set(unique(res.copies[end])) ⊆ Set([1, 2])
    # concentration homeostasis without replication: mRNA count scales with volume
    st_norep = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), threads = false)
    rn = simulate_population(tm, x0, 300, (0.0, 100.0); settings = st_norep, rng = Xoshiro(2))
    sn = final_snapshot(rn)
    v = sn.volume; c = sn.counts[:, 3]
    big = v .> 1.6; small = v .< 1.3
    @test abs(mean(c[big]) / mean(c[small]) - mean(v[big]) / mean(v[small])) < 0.25
    @test abs(mean(c[big] ./ v[big]) / mean(c[small] ./ v[small]) - 1) < 0.15
    # with replication, large (replicated) cells carry more transcripts than volume scaling alone predicts
    sn2 = final_snapshot(res)
    v2 = sn2.volume; c2 = sn2.counts[:, 3]
    @test mean(c2[v2 .> 1.6]) / mean(c2[v2 .< 1.3]) > mean(v2[v2 .> 1.6]) / mean(v2[v2 .< 1.3])
    # deterministic given seed, independent of threading
    st2 = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), threads = true)
    r1 = simulate_population(tm, x0, 50, (0.0, 20.0); settings = st2, rng = Xoshiro(3))
    r2 = simulate_population(tm, x0, 50, (0.0, 20.0); settings = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), threads = false), rng = Xoshiro(3))
    @test r1.counts[end] == r2.counts[end] && r1.volume[end] == r2.volume[end]
    # free growth doubles roughly every ln2/λ
    st3 = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), control = FreeGrowth(max_cells = 5000), threads = false)
    r3 = simulate_population(tm, x0, 100, (0.0, 30.0); settings = st3, rng = Xoshiro(4))
    @test 500 < r3.popsize[end] < 1100
    # logistic control saturates below K
    st4 = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), control = LogisticGrowth(400.0; max_cells = 5000), threads = false)
    r4 = simulate_population(tm, x0, 100, (0.0, 120.0); settings = st4, rng = Xoshiro(5))
    @test 250 < r4.popsize[end] < 550
    # adder and timer size control run and divide
    for sc in (Adder(1.0; cv = 0.05), AgeTimer(10.0; cv = 0.05))
        r = simulate_population(tm, x0, 50, (0.0, 40.0); settings = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = sc, threads = false), rng = Xoshiro(6))
        @test count(==(:divided), r.lineage.fate) > 50
    end
    # partitioning: binomial split conserves molecules, beta-binomial is over-dispersed
    tm4 = telegraph_model(k_on = 1.0, k_off = 1.0, k_tx = 40.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
    xa = zeros(Int, 4); xb = zeros(Int, 4)
    x = [1, 0, 1000, 500]
    BM.partition!(xa, xb, x, 0.5, tm4, BinomialPartition(), [1], Xoshiro(7))
    @test xa .+ xb == [2, 0, 1000, 500] || (xa[1:2] == x[1:2] && xb[1:2] == x[1:2] && xa[3:4] .+ xb[3:4] == x[3:4])
    rng = Xoshiro(8)
    bin = [BM.partition_count(1000, 0.5, BinomialPartition(), rng) for _ in 1:2000]
    bb = [BM.partition_count(1000, 0.5, BetaBinomialPartition(ρ = 0.05), rng) for _ in 1:2000]
    @test var(bb) > 3 * var(bin)
    # replicated promoter is split one copy per daughter
    xr = [1, 1, 10, 10]
    BM.partition!(xa, xb, xr, 0.5, tm4, BinomialPartition(), [1], Xoshiro(9))
    @test xa[1] + xa[2] == 1 && xb[1] + xb[2] == 1
    @test_throws ArgumentError simulate_population(tm, zeros(Int, 4), 10, (0.0, 1.0))
end
