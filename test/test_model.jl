@testset "model construction" begin
    m = birth_death_model(k = 3.0, γ = 0.5)
    @test nspecies(m) == 1 && nreactions(m) == 2 && nparams(m) == 2
    @test speciesnames(m) == [:X]
    @test stoichiometry(m) == [1 -1]
    @test m.p0[paramindex(m, :k)] == 3.0
    @test set_params(m; γ = 2.0)[paramindex(m, :γ)] == 2.0
    @test_throws ArgumentError set_params(m; nope = 1.0)
    @test initial_state(m; X = 4) == [4]

    tm = telegraph_model(k_on = 0.1, k_off = 0.2, k_tx = 5.0, k_dm = 1.0, k_tl = 2.0, k_dp = 0.1)
    @test nspecies(tm) == 4
    @test tm.is_promoter == [true, true, false, false]
    @test tm.promoter_groups == [[1, 2]]
    @test speciesindex(tm, :protein) == 4
    S = stoichiometry(tm)
    @test S[:, 1] == [-1, 1, 0, 0]          # activation
    @test S[:, 3] == [0, 0, 1, 0]           # transcription leaves G_on unchanged
    @test S[:, 5] == [0, 0, 0, 1]           # translation
    # dependency graph: firing activation must update activation, inactivation, transcription
    @test sort(tm.depgraph[1]) == [1, 2, 3]
    @test tm.depgraph[4] == [4, 5]           # mRNA degradation affects itself and translation
    @test tm.depgraph[3] == [4, 5]           # transcription affects mRNA degradation and translation

    # propensities
    x = initial_state(tm; G_on = 1, mRNA = 3, protein = 10)
    a = propensities(tm, x)
    @test a ≈ [0.0, 0.2, 5.0, 3.0, 6.0, 1.0]
    a2 = propensities(tm, x; V = 2.0)
    @test a2[3] ≈ 10.0                        # transcription ∝ V
    @test a2[4] ≈ 3.0                         # first order unchanged
    @test propensities(tm, x; copies = 2) == a  # promoter model: copies handled by promoter counts

    # second-order mass action and volume scaling
    dim = ReactionModel([Reaction("dimerise", [:A => 2], [:D], MassAction(:kd)),
                         Reaction("bind", [:A, :B], [:C], MassAction(:kb))]; params = (kd = 1.0, kb = 2.0))
    x = initial_state(dim; A = 5, B = 3)
    a = propensities(dim, x; V = 2.0)
    @test a[1] ≈ 1.0 * 10 / 2                 # binomial(5,2) / V
    @test a[2] ≈ 2.0 * 15 / 2
    @test dim.hor[speciesindex(dim, :A)] == 2 && dim.mrr[speciesindex(dim, :A)] == 2

    # Hill kinetics with literal and named parameters
    h = ReactionModel([Reaction("tx", [], [:M], Hill(:k; activators = [:A], inhibitors = [:R], K = [:KA, 4.0], n = 2.0, basal = 0.1)),
                       Reaction("deg", [:M], [], MassAction(:g))];
                      params = (k = 10.0, g = 1.0, KA = 2.0))
    x = initial_state(h; A = 2, R = 0)
    a = propensities(h, x; V = 1.0)
    @test a[1] ≈ 10.0 * (0.1 + 0.9 * (4 / (4 + 4)) * 1.0)
    x = initial_state(h; A = 2, R = 4)
    @test propensities(h, x)[1] ≈ 10.0 * (0.1 + 0.9 * 0.5 * (16 / (16 + 16)))
    hor = ReactionModel([Reaction("tx", [], [:M], Hill(:k; activators = [:A, :B], K = 1.0, n = 1.0, logic = :or)),
                         Reaction("deg", [:M], [], MassAction(:g))]; params = (k = 1.0, g = 1.0))
    @test propensities(hor, initial_state(hor; A = 1, B = 1))[1] ≈ 1 - 0.25
    # Custom kinetics
    c = ReactionModel([Reaction("in", [], [:X], Custom((x, V, θ, t) -> θ[1] * (1 + sin(t)); params = [:k])),
                       Reaction("out", [:X], [], MassAction(:g))]; params = (k = 2.0, g = 1.0))
    @test propensities(c, [0]; t = 0.0)[1] ≈ 2.0
    @test_throws ArgumentError ReactionModel([Reaction("in", [], [:X], MassAction(:missing))])
    @test_throws ArgumentError Reaction("bad", [], [:X], MassAction(:k); volume = :weird)
end
