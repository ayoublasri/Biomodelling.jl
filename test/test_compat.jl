@testset "v1 compatibility layer" begin
    k1 = 1.0; k2 = 0.2; k3 = 10.0; k4 = 0.1
    reaction1 = (name = "transcription", rate = k1, reactants = [:NULL], products = [:mRNA], coeff_rea = [1], coeff_pro = [1])
    reaction2 = (name = "mRNA decay", rate = k2, reactants = [:mRNA], products = [:NULL], coeff_rea = [1], coeff_pro = [1])
    reaction3 = (name = "translation", rate = k3, reactants = [:mRNA], products = [:mRNA, :protein], coeff_rea = [1], coeff_pro = [1, 1])
    reaction4 = (name = "protein decay", rate = k4, reactants = [:protein], products = [:NULL], coeff_rea = [1], coeff_pro = [1])
    model = (reaction1, reaction2, reaction3, reaction4)
    initiale_population = [:NULL 0; :mRNA 5; :protein 100]
    data = Donne(model, initiale_population, 2000.0, 0.5, 20, 0.03, 0.03)
    @test data.M == 4 && data.N == 3 && data.species == [:NULL, :mRNA, :protein] && data.growth_rate == 0.03
    T, X = ssa(data; rng = Xoshiro(1))
    @test size(X) == (4001, 3) && all(X[:, 1] .== 0)
    @test abs(mean(X[1000:end, 2]) - 5.0) < 0.4
    @test abs(mean(X[1000:end, 3]) - 500.0) < 50.0
    T2, X2 = tauleapswitch(data, 100; rng = Xoshiro(2))
    @test size(X2) == (4001, 3) && abs(mean(X2[1000:end, 3]) - 500.0) < 50.0
    T3, X3 = adaptive_tauleap(data; rng = Xoshiro(3))
    @test size(X3) == (4001, 3) && abs(mean(X3[1000:end, 3]) - 500.0) < 50.0
    data_short = Donne(model, initiale_population, 100.0, 0.5, 20, 0.03, 0.03)
    t, V, Xp = exponential_growth(data_short, 0.03, ssa, 1.0; rng = Xoshiro(4))
    @test size(V) == (201, 20) && size(Xp) == (201, 20, 3)
    @test all(1.0 .<= V[1, :] .<= 2.0)
    # second-order reaction from the original test suite
    Reaction1 = (name = "trans", rate = 4e-5, reactants = [:A; :B], products = [:A], coeff_rea = [2; 1], coeff_pro = [3])
    Reaction2 = (name = "birth", rate = 50, reactants = [:NULL], products = [:A], coeff_rea = [1], coeff_pro = [1])
    Reaction3 = (name = "death", rate = 10, reactants = [:A], products = [:NULL], coeff_rea = [1], coeff_pro = [1])
    Reaction4 = (name = "transf", rate = 25, reactants = [:NULL], products = [:B], coeff_rea = [1], coeff_pro = [1])
    d2 = Donne((Reaction1, Reaction2, Reaction3, Reaction4), [:NULL 0; :A 10; :B 10], 100, 1.0, 10, 0.38)
    S = stoichiometry(d2.rmodel)
    @test S[:, 1] == [1, -1]      # 2A + B → 3A: net A +1, B -1
    @test d2.rmodel.reactant_stoich[1] == [(1, 2), (2, 1)]
    # regulated reactions and growth-rate estimation
    act = (name = "activation", rate = [2.0, 2.0, 5.0], reactants = [:A], products = [:A, :C], coeff_rea = [1], coeff_pro = [1, 1])
    d3 = Donne((Reaction2, Reaction3, act, (name = "cdeath", rate = 1.0, reactants = [:C], products = [:NULL], coeff_rea = [1], coeff_pro = [1])),
               [:NULL 0; :A 5; :C 0], 50.0, 0.5, 5, 0.03)
    a = propensities(d3.rmodel, [5, 0])
    @test a[3] ≈ 2.0 * 25 / (25 + 25)
    g = growth_estimate(hcat(0:10, zeros(11), zeros(11), zeros(11), 100 .* exp.(0.2 .* (0:10))))
    @test g ≈ 0.2 atol = 1e-6
    @test_throws ErrorException Donne(model, zeros(2, 2), initiale_population, 10.0, 0.1, 5, 0.03)
    # the deprecated `MDR1` birth-death check from the original tests
    Rb = (name = "birth", rate = 100.0, reactants = [:NULL], products = [:MDR1], coeff_rea = [1], coeff_pro = [1])
    Rd = (name = "death", rate = 1.0, reactants = [:MDR1], products = [:NULL], coeff_rea = [1], coeff_pro = [1])
    d4 = Donne((Rb, Rd), [:NULL 0; :MDR1 100], 500.0, 1.0, 10, 0.38)
    C, D = tauleapswitch(d4, 100; rng = Xoshiro(5))
    @test abs(mean(D[:, 2]) - 100.0) < 2.0
end
