# Figure 5: ground-truth benchmarks — memory genes, GRN inference under division, imputation, perturbations.
include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra
const λ = log(2) / 20

# (a) memory genes: 40 independent telegraph genes with switching rates spanning three decades
ks = exp10.(range(-3, 0; length = 40))
rx = Reaction[]; params = Dict{Symbol,Float64}(); groups = Vector{Vector{Symbol}}()
for (j, k) in enumerate(ks)
    off = Symbol("G$(j)_off"); on = Symbol("G$(j)_on"); mr = Symbol("mRNA_$j")
    push!(rx, Reaction("act_$j", [off], [on], MassAction(Symbol("k_$j"))))
    push!(rx, Reaction("inact_$j", [on], [off], MassAction(Symbol("k_$j"))))
    push!(rx, Reaction("tx_$j", [on], [on, mr], MassAction(:k_tx); volume = :proportional))
    push!(rx, Reaction("deg_$j", [mr], [], MassAction(:k_dm)))
    params[Symbol("k_$j")] = k
    push!(groups, [off, on])
end
params[:k_tx] = 20.0; params[:k_dm] = 0.5
mm = ReactionModel(rx; params = params, promoters = groups)
x0 = zeros(Int, nspecies(mm)); for g in mm.promoter_groups; x0[g[1]] = 1; end
st = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), control = FreeGrowth(max_cells = 20000), record_every = 1000)
r = simulate_population(mm, x0, 40, (0.0, 120.0); settings = st, rng = Xoshiro(1))       # 40 founder clones, ~6 doublings
sn = final_snapshot(r)
mi = [speciesindex(mm, Symbol("mRNA_$j")) for j in 1:40]
true_scores = clonal_variance_scores(sn.counts[:, mi] ./ sn.volume, sn.clone; n_perm = 300, rng = Xoshiro(2))
seq = sequence(sn.counts[:, mi], SeqProtocol(capture = 0.2, capture_cv = 0.3); rng = Xoshiro(3))
Ynorm = seq.Y ./ max.(sum(seq.Y; dims = 2), 1) .* median(sum(seq.Y; dims = 2))
seq_scores = clonal_variance_scores(Ynorm, sn.clone; n_perm = 300, rng = Xoshiro(4))
save_csv("fig5a_memory_genes.csv", ["gene", "k_switch", "memory_time", "memory_generations", "score_true", "p_true", "score_seq", "p_seq", "mean_count"],
         hcat(1:40, ks, 1 ./ (2ks), 1 ./ (2ks) ./ 20, true_scores.score, true_scores.pvalue, seq_scores.score, seq_scores.pvalue, vec(mean(sn.counts[:, mi]; dims = 1))))
println("(a) cells = ", size(sn.counts, 1), " clones = ", length(unique(sn.clone)))

# (b) GRN inference: same network at fixed volume, in a dividing population (counts / concentrations), and after sequencing
grn, adj = random_grn(30; n_activations = 30, n_inhibitions = 15, telegraph = false, k_tx = (5.0, 30.0), K = (2.0, 10.0), n = 2.0, basal = 0.05, rng = Xoshiro(5))
save_csv("fig5b_adjacency.csv", ["gene_$j" for j in 1:30], adj)
xg = zeros(Int, nspecies(grn))
datasets = Dict{String,Matrix{Float64}}()
E = ensemble_final(grn, xg, 80.0, 1500; rng = Xoshiro(6))
datasets["fixed_volume"] = float.(E)
stp = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), replication = Replication(0.5), record_every = 1000)
rp = simulate_population(grn, xg, 1500, (0.0, 100.0); settings = stp, rng = Xoshiro(7))
snp = final_snapshot(rp)
datasets["population_counts"] = float.(snp.counts)
datasets["population_concentration"] = snp.counts ./ snp.volume
sq = sequence(snp.counts, SeqProtocol(capture = 0.15, capture_cv = 0.3); rng = Xoshiro(8))
datasets["sequenced_counts"] = float.(sq.Y)
ls = vec(sum(sq.Y; dims = 2))
datasets["sequenced_normalized"] = sq.Y ./ max.(ls, 1) .* median(ls)
for (name, Y) in datasets
    save_csv("fig5b_data_$name.csv", ["gene_$j" for j in 1:30], Y)
end
save_csv("fig5b_population_metadata.csv", ["volume", "age", "copies", "library_size"], hcat(snp.volume, snp.age, snp.copies, ls))
# correlation-based and PIDC scores
function spearman_matrix(Y)
    R = mapslices(x -> StatsBase.tiedrank(x), Y; dims = 1)
    cor(R)
end
for (name, Y) in datasets
    Z = log1p.(Y)
    C = abs.(cor(Z)); C[diagind(C)] .= 0
    S = abs.(spearman_matrix(Y)); S[diagind(S)] .= 0
    save_csv("fig5b_scores_pearson_$name.csv", ["gene_$j" for j in 1:30], C)
    save_csv("fig5b_scores_spearman_$name.csv", ["gene_$j" for j in 1:30], S)
end
# PIDC scores are produced by pidc_scores.jl (separate environment, see paper/README.md)

# (d) perturbation ground truth: knockdowns of the five strongest regulators
g20, adj20 = random_grn(20; n_activations = 22, n_inhibitions = 10, telegraph = false, k_tx = (5.0, 30.0), K = (2.0, 10.0), n = 2.0, basal = 0.05, rng = Xoshiro(9))
save_csv("fig5d_adjacency.csv", ["gene_$j" for j in 1:20], adj20)
x20 = zeros(Int, nspecies(g20))
stc = PopulationSettings(dt = 0.1, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), record_every = 1000)
ctrl = simulate_population(g20, x20, 800, (0.0, 100.0); settings = stc, rng = Xoshiro(10))
Yc = final_snapshot(ctrl).counts ./ final_snapshot(ctrl).volume
mu_c = vec(mean(Yc; dims = 1))
Rc = cor(log1p.(Yc))
regs = sortperm(vec(sum(adj20 .!= 0; dims = 2)); rev = true)[1:5]
rows = Any[]
for j in regs
    pert = Perturbation(; gene_perturbations = [GenePerturbation(Symbol("k_tx_$j"), 0.05; fraction = 1.0, t_start = 0.0)])
    kd = simulate_population(g20, x20, 800, (0.0, 100.0); settings = stc, perturbation = pert, rng = Xoshiro(10))
    Yk = final_snapshot(kd).counts ./ final_snapshot(kd).volume
    mu_k = vec(mean(Yk; dims = 1))
    for i in 1:20
        push!(rows, [j, i, mu_c[i], mu_k[i], log2((mu_k[i] + 0.1) / (mu_c[i] + 0.1)), Rc[i, j], adj20[j, i]])
    end
end
save_csv("fig5d_knockdowns.csv", ["knocked_gene", "gene", "mean_control", "mean_knockdown", "log2fc_true", "correlation_with_knocked", "direct_edge"], permutedims(reduce(hcat, rows)))
println("fig5 (julia part) done; run benchmarks_python.py for GENIE3, imputation and metrics")
