# Figure 3: emergent memory, size scaling and cell-cycle effects in growing and dividing cells.
include(joinpath(@__DIR__, "common.jl"))

λ = log(2) / 20                                    # 20 time-unit doubling time
base_settings(; kwargs...) = PopulationSettings(; dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05),
                                                partitioning = BinomialPartition(σ = 0.02), kwargs...)
mk(k) = telegraph_model(k_on = k, k_off = k, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2)

# (a) single-lineage traces of volume, mRNA and protein for a slowly switching gene
slow = mk(0.005)
x0 = initial_state(slow; G_off = 1)
res = simulate_population(slow, x0, 300, (0.0, 240.0); settings = base_settings(replication = Replication(0.5)), rng = Xoshiro(1))
# pick the surviving cell whose ancestral line switched promoter state most often (an informative trace)
best = (-1, res.ids[end][1])
for id in res.ids[end][1:min(60, end)]
    t = follow_lineage(res, id)
    nsw = count(i -> (t.X[i, 2] > 0) != (t.X[i - 1, 2] > 0), 2:size(t.X, 1))
    nsw > best[1] && (best = (nsw, id))
end
tr = follow_lineage(res, best[2])
save_csv("fig3a_trace.csv", ["t", "V", "G_on", "mRNA", "protein", "cell_id"], hcat(tr.t, tr.V, tr.X[:, 2], tr.X[:, 3], tr.X[:, 4], tr.ids))

# (b) count-volume scaling with and without gene replication (fast gene → snapshot at steady state)
fast = mk(1.0)
for (tag, rep) in (("no_replication", nothing), ("replication", Replication(0.5)))
    r = simulate_population(fast, initial_state(fast; G_off = 1), 2000, (0.0, 120.0); settings = base_settings(replication = rep, record_every = 1200), rng = Xoshiro(2))
    sn = final_snapshot(r)
    save_csv("fig3b_scaling_$tag.csv", ["volume", "mRNA", "protein", "copies", "age"], hcat(sn.volume, sn.counts[:, 3], sn.counts[:, 4], sn.copies, sn.age))
end

# (c) heritability versus promoter switching rate (mother-daughter, sisters, cousins)
rows = Any[]
for k in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 1.0)
    m = mk(k)
    r = simulate_population(m, initial_state(m; G_off = 1), 400, (0.0, 240.0); settings = base_settings(record_every = 2400), rng = Xoshiro(3))
    md = heritability(r, :protein; relation = :mother_daughter); ss = heritability(r, :protein; relation = :sisters); cc = heritability(r, :protein; relation = :cousins)
    mdm = heritability(r, :mRNA; relation = :mother_daughter)
    push!(rows, [k, md.r, ss.r, cc.r, mdm.r, md.n])
    @printf("k=%.3f  mother-daughter r=%.3f sisters r=%.3f cousins r=%.3f (n=%d)\n", k, md.r, ss.r, cc.r, md.n)
end
save_csv("fig3c_heritability.csv", ["k_switch", "mother_daughter", "sisters", "cousins", "mother_daughter_mRNA", "n_pairs"], permutedims(reduce(hcat, rows)))

# (d) memory timescale (time units) versus cell-cycle time for fixed switching rates
rows = Any[]
for Td in (10.0, 20.0, 40.0), k in (0.005, 0.02)
    m = mk(k)
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(log(2) / Td), size_control = Sizer(2.0; cv = 0.05), partitioning = BinomialPartition(σ = 0.02), record_every = 100_000)
    r = simulate_population(m, initial_state(m; G_off = 1), 400, (0.0, 12Td); settings = st, rng = Xoshiro(4))
    mt = memory_timescale(r, :protein; max_generations = 6)
    push!(rows, [Td, k, mt.generations, mt.time, mt.cycle_time, mt.correlations...])
end
save_csv("fig3d_memory_timescale.csv", vcat(["cycle_time_nominal", "k_switch", "tau_generations", "tau_time", "cycle_time_measured"], ["r_gen$g" for g in 1:6]), permutedims(reduce(hcat, rows)))

# (e) lineage versus population noise as partitioning noise and cycle-time noise vary
rows = Any[]
const_gene = telegraph_model(k_on = 5.0, k_off = 5.0, k_tx = 20.0, k_dm = 1.0, k_tl = 2.0, k_dp = 0.05)   # stable protein, diluted by growth
for σp in (0.0, 0.1), cvT in (0.0, 0.15, 0.3), part in ("binomial", "betabinomial")
    pt = part == "binomial" ? BinomialPartition(σ = σp) : BetaBinomialPartition(σ = σp, ρ = 0.1)
    st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = cvT), partitioning = pt)
    r = simulate_population(const_gene, initial_state(const_gene; G_off = 1), 400, (0.0, 200.0); settings = st, rng = Xoshiro(5))
    nd = noise_decomposition(r, :protein; n_lineages = 60, burnin = 0.3)
    push!(rows, [σp, cvT, part, nd.population, nd.lineage, nd.n_lineages])
end
save_csv("fig3e_noise.csv", ["partition_sigma", "cycle_cv", "partitioning", "cv2_population", "cv2_lineage", "n_lineages"], permutedims(reduce(hcat, rows)))

# (f) inferred burst parameters before and after gene replication (copy number 1 vs 2), and when the cycle is ignored
g = telegraph_model(k_on = 0.5, k_off = 1.5, k_tx = 40.0, k_dm = 1.0)
r = simulate_population(g, initial_state(g; G_off = 1), 3000, (0.0, 120.0); settings = base_settings(replication = Replication(0.5), record_every = 1200), rng = Xoshiro(6))
sn = final_snapshot(r)
mi = speciesindex(g, :mRNA)
rows = Any[]
for (tag, sel) in (("copies1", sn.copies .== 1), ("copies2", sn.copies .== 2), ("pooled", trues(length(sn.copies))))
    f = fit_telegraph(sn.counts[sel, mi])
    push!(rows, [tag, count(sel), f.k_on, f.k_off, f.k_tx, mean(sn.counts[sel, mi]), mean(sn.volume[sel])])
    @printf("%-8s n=%5d k_on=%.3f k_off=%.3f k_tx=%.2f\n", tag, count(sel), f.k_on, f.k_off, f.k_tx)
end
save_csv("fig3f_copynumber.csv", ["subset", "n_cells", "k_on", "k_off", "k_tx", "mean_mRNA", "mean_volume"], permutedims(reduce(hcat, rows)))
println("fig3 done")
