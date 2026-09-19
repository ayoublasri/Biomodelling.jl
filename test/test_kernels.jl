@testset "kernels: exactness and agreement" begin
    rng = Xoshiro(11)
    # birth-death → Poisson(k/γ)
    bd = birth_death_model(k = 12.0, γ = 1.0)
    pmf = pdf.(Poisson(12.0), 0:60)
    for (kern, tol) in ((DirectSSA(), 0.015), (HybridSSATau(0.02), 0.02), (AdaptiveTauLeap(), 0.02), (TauLeap(0.005), 0.02))
        E = ensemble_final(bd, [12], 15.0, 8000; kernel = kern, rng = Xoshiro(1))
        @test minimum(E) >= 0
        @test abs(mean(E[:, 1]) - 12.0) < 0.2
        @test ks_discrete(E[:, 1], pmf) < tol
    end
    # telegraph → Beta-Poisson (Peccoud–Ycart)
    tm = telegraph_model(k_on = 0.4, k_off = 0.6, k_tx = 12.0, k_dm = 1.0; volume_scaled = false)
    x0 = initial_state(tm; G_off = 1)
    E = ensemble_final(tm, x0, 30.0, 12000; rng = Xoshiro(2))
    c = E[:, speciesindex(tm, :mRNA)]
    pmf = telegraph_pmf(0:80, 0.4, 0.6, 12.0, 1.0)
    @test sum(pmf) ≈ 1.0 atol = 1e-8
    @test abs(mean(c) - 12.0 * 0.4) < 0.15
    @test ks_discrete(c, pmf) < 0.015
    @test abs(mean(E[:, speciesindex(tm, :G_on)]) - 0.4) < 0.02
    # bursty protein (two-stage with short-lived mRNA) → negative binomial
    a, b, γ = 1.5, 6.0, 1.0
    bp = bursty_protein_model(a = a, b = b, γ = γ)
    E = ensemble_final(bp, initial_state(bp), 25.0, 8000; rng = Xoshiro(3))
    prot = E[:, speciesindex(bp, :protein)]
    nb = NegativeBinomial(a / γ, 1 / (1 + b))
    @test abs(mean(prot) - a * b / γ) < 0.4
    @test ks_discrete(prot, pdf.(nb, 0:200)) < 0.03
    # single-cell trajectory recording and determinism
    tr = simulate(bd, [0], (0.0, 10.0); saveat = 0:1.0:10, rng = Xoshiro(4))
    tr2 = simulate(bd, [0], (0.0, 10.0); saveat = 0:1.0:10, rng = Xoshiro(4))
    @test tr.X == tr2.X && length(tr.t) == 11 && tr[:X] == tr.X[:, 1]
    E1 = ensemble_final(bd, [0], 5.0, 300; rng = Xoshiro(5), threads = false)
    E2 = ensemble_final(bd, [0], 5.0, 300; rng = Xoshiro(5), threads = true)
    @test E1 == E2
    # strict tau-leap errors instead of going negative
    @test_throws ErrorException ensemble_final(birth_death_model(k = 1.0, γ = 5.0), [1], 5.0, 400; kernel = TauLeap(0.5), rng = Xoshiro(6), threads = false)
    # adaptive tau-leap on a stiff dimerisation system stays non-negative and matches SSA mean
    dim = ReactionModel([Reaction("in", [], [:A], MassAction(:k)), Reaction("dimer", [:A => 2], [:D], MassAction(:kd)),
                         Reaction("undimer", [:D], [:A => 2], MassAction(:ku)), Reaction("out", [:D], [], MassAction(:g))];
                        params = (k = 50.0, kd = 0.05, ku = 1.0, g = 0.5))
    Es = ensemble_final(dim, [0, 0], 30.0, 2000; rng = Xoshiro(7))
    Ea = ensemble_final(dim, [0, 0], 30.0, 2000; kernel = AdaptiveTauLeap(ε = 0.02), rng = Xoshiro(8))
    @test minimum(Ea) >= 0
    @test abs(mean(Ea[:, 2]) - mean(Es[:, 2])) / mean(Es[:, 2]) < 0.05
end
