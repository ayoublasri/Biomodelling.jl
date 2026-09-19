@testset "inference" begin
    # telegraph pmf: normalisation, moments and limits
    p = telegraph_pmf(0:200, 0.5, 1.5, 20.0, 1.0)
    @test sum(p) ≈ 1.0 atol = 1e-8
    m = sum(p .* (0:200))
    @test m ≈ 20.0 * 0.25 atol = 1e-6
    v = sum(p .* (0:200) .^ 2) - m^2
    @test v ≈ m + 20.0^2 * 0.25 * 0.75 / (0.5 + 1.5 + 1.0) atol = 1e-5
    # constitutive limit: k_on ≫ k_off → Poisson
    pc = telegraph_pmf(0:60, 500.0, 0.5, 10.0, 1.0)
    @test maximum(abs.(pc .- pdf.(Poisson(10.0 * 500 / 500.5), 0:60))) < 1e-3
    @test telegraph_pmf(3, 0.5, 1.5, 20.0, 2.0) ≈ telegraph_pmf(3, 0.25, 0.75, 10.0, 1.0)
    # MLE recovers parameters from simulated data
    tm = telegraph_model(k_on = 0.4, k_off = 0.8, k_tx = 25.0, k_dm = 1.0; volume_scaled = false)
    E = ensemble_final(tm, initial_state(tm; G_off = 1), 30.0, 6000; rng = Xoshiro(1))
    c = E[:, speciesindex(tm, :mRNA)]
    fit = fit_telegraph(c)
    @test fit.converged
    @test abs(log(fit.k_tx / 25.0)) < 0.25
    @test abs(log(fit.k_on / 0.4)) < 0.5
    @test abs(log(fit.k_off / 0.8)) < 0.5
    @test telegraph_loglik(c, fit.k_on, fit.k_off, fit.k_tx) >= telegraph_loglik(c, 0.4, 0.8, 25.0) - 3.0
    # summaries and distances
    s = moment_summaries(c)
    @test length(s) == 4 && s[1] ≈ mean(c) && s[4] ≈ mean(c .== 0)
    @test summary_distance(s, s) == 0.0
    @test BM.wasserstein1(c, c) == 0.0 && BM.wasserstein1(c, c .+ 2) ≈ 2.0
    # ABC-SMC recovers a birth rate
    bd = birth_death_model(k = 20.0, γ = 1.0)
    obs = ensemble_final(bd, [0], 15.0, 2000; rng = Xoshiro(2))[:, 1]
    sobs = moment_summaries(obs)
    dist(θ, rng) = summary_distance(moment_summaries(ensemble_final(bd, [0], 15.0, 200; p = set_params(bd; k = θ[1]), rng = rng, threads = false)[:, 1]), sobs)
    res = abc_smc(dist, [Uniform(1.0, 80.0)]; n_particles = 80, generations = 5, rng = Xoshiro(3), names = [:k])
    pm = posterior_mean(res)[1]
    @test abs(pm - 20.0) < 3.0
    lo, hi = credible_interval(res, :k)
    @test lo < 20.0 < hi && hi - lo < 20.0
    @test length(res.epsilons) == 5 && issorted(res.epsilons[2:end]; rev = true)
end
