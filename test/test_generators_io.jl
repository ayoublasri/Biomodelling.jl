using HDF5
@testset "generators and io" begin
    m, adj = random_grn(6; n_activations = 5, n_inhibitions = 3, rng = Xoshiro(1))
    @test count(==(1), adj) == 5 && count(==(-1), adj) == 3
    @test all(adj[i, i] == 0 for i in 1:6)
    @test nspecies(m) == 6 && nreactions(m) == 12
    mt, adjt = random_grn(5; n_activations = 4, n_inhibitions = 2, telegraph = true, topology = :scale_free, k_tx = (5.0, 20.0), rng = Xoshiro(2))
    @test nspecies(mt) == 15 && length(mt.promoter_groups) == 5
    x0 = zeros(Int, nspecies(mt)); for g in mt.promoter_groups; x0[g[1]] = 1; end
    res = simulate_population(mt, x0, 60, (0.0, 20.0); settings = PopulationSettings(dt = 0.2, threads = false), rng = Xoshiro(3))
    @test size(res.counts[end]) == (60, 15)
    @test_throws ArgumentError random_grn(3; n_activations = 5, n_inhibitions = 5)
    # csv and newick
    dir = mktempdir()
    sn = final_snapshot(res)
    f1 = write_counts_csv(joinpath(dir, "c.csv"), sn.counts, sn.species, sn.ids)
    @test countlines(f1) == 61
    f2 = write_metadata_csv(joinpath(dir, "m.csv"), sn)
    @test countlines(f2) == 61
    f3 = write_lineage_csv(joinpath(dir, "l.csv"), res.lineage)
    @test countlines(f3) == length(res.lineage) + 1
    nw = newick(res.lineage; t_end = 20.0)
    @test endswith(nw, ";") && occursin("c1:", nw)
    # h5ad through the HDF5 extension
    path = joinpath(dir, "d.h5ad")
    write_h5ad(path, sn.counts; genes = sn.species, cells = sn.ids,
               obs = (volume = sn.volume, clone = sn.clone, perturbed = sn.perturbed, label = fill("x", 60)),
               var = (kind = String.(sn.species),), obsm = (X_true = Float64.(sn.counts),), uns = (seed = 3, note = "test"))
    h5open(path) do f
        @test Set(keys(f)) ⊇ Set(["X", "obs", "var", "obsm", "uns"])
        @test size(f["X"]) == (15, 60)          # row-major (n_obs, n_vars) for anndata
        @test attrs(f["obs"])["encoding-type"] == "dataframe"
        @test Set(attrs(f["obs"])["column-order"]) == Set(["volume", "clone", "perturbed", "label"])
        @test length(read(f["obs"]["_index"])) == 60
        @test read(f["uns"]["note"]) == "test"
    end
end
