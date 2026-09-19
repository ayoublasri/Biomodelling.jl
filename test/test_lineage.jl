@testset "lineage statistics" begin
    λ = log(2) / 20
    slow = telegraph_model(k_on = 0.003, k_off = 0.003, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
    fast = telegraph_model(k_on = 3.0, k_off = 3.0, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)
    x0 = initial_state(slow; G_off = 1)
    st = PopulationSettings(dt = 0.2, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), threads = false)
    rs = simulate_population(slow, x0, 300, (0.0, 200.0); settings = st, rng = Xoshiro(1))
    rf = simulate_population(fast, x0, 300, (0.0, 200.0); settings = st, rng = Xoshiro(1))
    hs = heritability(rs, :protein; relation = :mother_daughter)
    hf = heritability(rf, :protein; relation = :mother_daughter)
    @test hs.n > 100 && hf.n > 100
    @test hs.r > 0.6                 # slow promoter switching → strong memory
    @test hs.r > hf.r + 0.3          # fast switching → little memory
    ss = heritability(rs, :protein; relation = :sisters)
    @test ss.r > 0.5
    ac = lineage_autocorrelation(rs, :protein; max_generations = 4)
    @test ac[1] > ac[3] > -0.2
    mt = memory_timescale(rs, :protein; max_generations = 4)
    @test mt.generations > 1.0
    @test abs(mt.cycle_time - 20.0) < 2.0
    # lineage traces run from the start to the end of the simulation
    tr = follow_lineage(rs, rs.ids[end][1])
    @test tr.t[1] == rs.t[1] && tr.t[end] == rs.t[end] && length(unique(tr.ids)) >= 5
    nd = noise_decomposition(rs, :protein; n_lineages = 10)
    @test isfinite(nd.population) && isfinite(nd.lineage)
    # tree utilities
    lt = rs.lineage
    sp = sister_pairs(lt); md = mother_daughter_pairs(lt)
    @test all(lt.parent[a] == lt.parent[b] for (a, b) in sp)
    @test all(lt.parent[d] == m for (m, d) in md)
    @test lineage_of(lt, sp[1][1])[end] in lt.id[lt.parent .== 0]
    nw = newick(lt; founder = 1, t_end = 200.0)
    @test startswith(nw, "(") || startswith(nw, "c1:")
    @test endswith(nw, ";")
    # clonal variance scores separate memory from non-memory genes
    cvs = clonal_variance_scores(rs; n_perm = 100, rng = Xoshiro(2))
    cvf = clonal_variance_scores(rf; n_perm = 100, rng = Xoshiro(2))
    pi = speciesindex(slow, :protein)
    @test cvs.score[pi] > 2.0
    @test cvs.score[pi] > cvf.score[pi]
    @test cvs.pvalue[pi] < 0.05
    # fluctuation test: memory gives super-binomial clone-to-clone variance
    ft = fluctuation_test(slow, x0, sn -> mean(sn.counts[:, pi] ./ sn.volume .> 50.0); n_clones = 20, generations = 4, settings = st, rng = Xoshiro(3))
    ftf = fluctuation_test(fast, x0, sn -> mean(sn.counts[:, pi] ./ sn.volume .> 50.0); n_clones = 20, generations = 4, settings = st, rng = Xoshiro(3))
    @test ft.ratio > 3.0
    @test ft.ratio > ftf.ratio
    λs = covariance_eigenspectrum(rs.counts[end])
    @test issorted(λs; rev = true) && length(λs) <= 4
    @test powerlaw_tail_exponent(λs; k = 3) < 0
end

@testset "kin pairs" begin
    lt = LineageTable()
    function add!(lt, id, parent, gen)
        push!(lt.id, id); push!(lt.parent, parent); push!(lt.clone, 1); push!(lt.generation, gen)
        push!(lt.birth_time, float(gen)); push!(lt.end_time, NaN); push!(lt.fate, :alive)
        push!(lt.V_birth, 1.0); push!(lt.V_end, NaN); push!(lt.x_birth, Int[]); push!(lt.x_end, Int[]); push!(lt.perturbed, false)
    end
    # founder 1 -> 2,3 ; 2 -> 4,5 ; 3 -> 6,7 ; 4 -> 8,9 ; 6 -> 10,11
    add!(lt, 1, 0, 0)
    for (id, p, gen) in ((2, 1, 1), (3, 1, 1), (4, 2, 2), (5, 2, 2), (6, 3, 2), (7, 3, 2), (8, 4, 3), (9, 4, 3), (10, 6, 3), (11, 6, 3))
        add!(lt, id, p, gen)
    end
    @test Set(kin_pairs(lt, 1)) == Set([(2, 3), (4, 5), (6, 7), (8, 9), (10, 11)])
    @test Set(kin_pairs(lt, 2)) == Set([(4, 6), (4, 7), (5, 6), (5, 7)])
    @test Set(kin_pairs(lt, 3)) == Set([(8, 10), (8, 11), (9, 10), (9, 11)])
    @test kin_pairs(lt, 2) == cousin_pairs(lt)
end
