@testset "observation models" begin
    rng = Xoshiro(1)
    counts = rand(rng, 0:200, 500, 6)
    r = sequence(counts, SeqProtocol(capture = 0.3, capture_cv = 0.2); rng = Xoshiro(2))
    @test size(r.Y) == size(counts) && all(r.Y .<= counts)
    @test abs(mean(r.Y) / mean(counts) - 0.3) < 0.03
    @test abs(mean(r.capture) - 0.3) < 0.02
    r2 = sequence(counts, SeqProtocol(capture = 0.3, capture_cv = 0.0, depth = 50.0); rng = Xoshiro(3))
    @test abs(mean(sum(r2.Y; dims = 2)) - 50.0) < 3.0
    r3 = sequence(counts, SeqProtocol(capture = 1.0, capture_cv = 0.0, dropout = 0.5); rng = Xoshiro(4))
    @test abs(mean(r3.Y .== 0) - 0.5) < 0.05
    r4 = sequence(counts, SeqProtocol(capture = 0.3, capture_cv = 0.0, batch_cv = 0.5, n_batches = 3); rng = Xoshiro(5))
    @test length(unique(r4.batch)) == 3
    @test_throws ArgumentError sequence(counts, SeqProtocol(capture = 0.5, capture_cv = 2.0))
    f = smfish(counts; efficiency = 0.9, rng = Xoshiro(6))
    @test abs(mean(f) / mean(counts) - 0.9) < 0.02
    # sampling and time-lapse from a population result
    tm = telegraph_model(k_on = 0.2, k_off = 0.2, k_tx = 20.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
    res = simulate_population(tm, initial_state(tm; G_off = 1), 100, (0.0, 30.0);
                              settings = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / 20), threads = false), rng = Xoshiro(7))
    sn = sample_cells(res; n = 40, rng = Xoshiro(8))
    @test size(sn.counts, 1) == 40 && length(sn.volume) == 40 && length(sn.clone) == 40
    tl = timelapse(res, :protein; interval = 1.0, σ = 0.5, n_lineages = 5, rng = Xoshiro(9))
    @test length(tl) == 5 && length(tl[1].t) == 31
    @test all(t -> t.t[1] == 0.0 && isapprox(t.t[end], 30.0; atol = 1e-6), tl)
    tl2 = timelapse(res, :protein; interval = 2.0, σ = -0.1, n_lineages = 2, normalize = :count)
    @test length(tl2[1].t) == 16
    @test BM.library_size(sn.counts) == vec(sum(sn.counts; dims = 2))
end
