# Figure 6: inference — ABC-SMC versus the exact likelihood, the bias of ignoring division, drug parameter recovery.
include(joinpath(@__DIR__, "common.jl"))
const λ = log(2) / 20
truth = (k_on = 0.3, k_off = 0.6, k_tx = 20.0)

# (a) snapshot of non-dividing cells: exact MLE and ABC-SMC agree
tm = telegraph_model(; truth..., k_dm = 1.0, volume_scaled = false)
x0 = initial_state(tm; G_off = 1)
mi = speciesindex(tm, :mRNA)
data = ensemble_final(tm, x0, 40.0, 2000; rng = Xoshiro(1))[:, mi]
mle = fit_telegraph(data)
sobs = moment_summaries(data)
dist_a(θ, rng) = summary_distance(moment_summaries(ensemble_final(tm, x0, 40.0, 400; p = set_params(tm; k_on = θ[1], k_off = θ[2], k_tx = θ[3]), rng = rng, threads = false)[:, mi]), sobs)
prior = [LogUniform(0.01, 10.0), LogUniform(0.01, 10.0), LogUniform(1.0, 200.0)]
abc_a = abc_smc(dist_a, prior; n_particles = 200, generations = 7, rng = Xoshiro(2), names = [:k_on, :k_off, :k_tx], verbose = true)
save_csv("fig6a_particles.csv", ["k_on", "k_off", "k_tx", "weight"], hcat(abc_a.particles, abc_a.weights))
save_csv("fig6a_schedule.csv", ["generation", "epsilon", "acceptance"], hcat(1:length(abc_a.epsilons), abc_a.epsilons, abc_a.acceptance))
save_kv("fig6a_summary.csv", ["true_k_on" => truth.k_on, "true_k_off" => truth.k_off, "true_k_tx" => truth.k_tx,
        "mle_k_on" => mle.k_on, "mle_k_off" => mle.k_off, "mle_k_tx" => mle.k_tx,
        "abc_mean_k_on" => posterior_mean(abc_a)[1], "abc_mean_k_off" => posterior_mean(abc_a)[2], "abc_mean_k_tx" => posterior_mean(abc_a)[3],
        "abc_ci_k_tx_lo" => credible_interval(abc_a, :k_tx)[1], "abc_ci_k_tx_hi" => credible_interval(abc_a, :k_tx)[2], "final_epsilon" => abc_a.epsilons[end]])
@printf("(a) MLE k_on %.3f k_off %.3f k_tx %.2f | ABC %.3f %.3f %.2f\n", mle.k_on, mle.k_off, mle.k_tx, posterior_mean(abc_a)...)

# (b) snapshot of a growing, dividing population: the naive telegraph fit is biased, the division-aware ABC is not
tmv = telegraph_model(; truth..., k_dm = 1.0, volume_scaled = true)
st = PopulationSettings(dt = 0.2, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), partitioning = BinomialPartition(σ = 0.02),
                        replication = Replication(0.5), record_every = 10_000, track_lineage = false)
pop = simulate_population(tmv, x0, 2000, (0.0, 100.0); settings = st, rng = Xoshiro(3))
snp = final_snapshot(pop)
counts_pop = snp.counts[:, mi]
naive = fit_telegraph(counts_pop)
sobs_b = moment_summaries(counts_pop)
function dist_b(θ, rng)
    p = set_params(tmv; k_on = θ[1], k_off = θ[2], k_tx = θ[3])
    r = simulate_population(tmv, x0, 200, (0.0, 80.0); settings = st, p = p, rng = rng)
    summary_distance(moment_summaries(final_snapshot(r).counts[:, mi]), sobs_b)
end
abc_b = abc_smc(dist_b, prior; n_particles = 120, generations = 5, rng = Xoshiro(4), names = [:k_on, :k_off, :k_tx], verbose = true)
save_csv("fig6b_particles.csv", ["k_on", "k_off", "k_tx", "weight"], hcat(abc_b.particles, abc_b.weights))
save_csv("fig6b_schedule.csv", ["generation", "epsilon", "acceptance"], hcat(1:length(abc_b.epsilons), abc_b.epsilons, abc_b.acceptance))
save_kv("fig6b_summary.csv", ["true_k_on" => truth.k_on, "true_k_off" => truth.k_off, "true_k_tx" => truth.k_tx,
        "naive_k_on" => naive.k_on, "naive_k_off" => naive.k_off, "naive_k_tx" => naive.k_tx,
        "abc_mean_k_on" => posterior_mean(abc_b)[1], "abc_mean_k_off" => posterior_mean(abc_b)[2], "abc_mean_k_tx" => posterior_mean(abc_b)[3],
        "abc_ci_k_tx_lo" => credible_interval(abc_b, :k_tx)[1], "abc_ci_k_tx_hi" => credible_interval(abc_b, :k_tx)[2],
        "mean_volume" => mean(snp.volume), "mean_copies" => mean(snp.copies), "final_epsilon" => abc_b.epsilons[end]])
@printf("(b) naive k_on %.3f k_off %.3f k_tx %.2f | ABC %.3f %.3f %.2f\n", naive.k_on, naive.k_off, naive.k_tx, posterior_mean(abc_b)...)

# (c) drug parameters from a kill curve and sister fate concordance
rm_ = telegraph_model(k_on = 0.01, k_off = 0.01, k_tx = 30.0, k_dm = 1.0, k_tl = 4.0, k_dp = 0.2, gene = :R, mrna = :mRNA, protein = :P)
xr = initial_state(rm_; R_off = 1)
truth_c = (h_max = 0.5, K = 150.0)
mkpert(h, K) = Perturbation(PiecewiseDose([0.0, 40.0], [0.0, 1.0]); effects = [DeathHazard(h_max = h, EC50 = 0.5, m = 2.0, protect = :P, K = K, q = 4.0)])
stc = PopulationSettings(dt = 0.2, growth = ExponentialGrowth(λ), size_control = Sizer(2.0; cv = 0.05), control = FreeGrowth(max_cells = 8000), record_every = 25)
function summaries_c(r)
    i0 = argmin(abs.(r.t .- 40.0))
    surv = log.(max.(r.popsize[i0:end], 1.0) ./ r.popsize[i0])
    lt = r.lineage
    same = 0; n = 0
    for (a, b) in sister_pairs(lt)
        (20.0 <= lt.birth_time[a] < 40.0 && 20.0 <= lt.birth_time[b] < 40.0) || continue
        n += 1; same += (lt.fate[a] == :died) == (lt.fate[b] == :died)
    end
    vcat(surv, [n == 0 ? 0.5 : same / n])
end
obs_c = simulate_population(rm_, xr, 400, (0.0, 100.0); settings = stc, perturbation = mkpert(truth_c...), rng = Xoshiro(5))
sobs_c = summaries_c(obs_c)
function dist_c(θ, rng)
    r = simulate_population(rm_, xr, 200, (0.0, 100.0); settings = stc, perturbation = mkpert(θ[1], θ[2]), rng = rng)
    s = summaries_c(r)
    sqrt(mean((s .- sobs_c) .^ 2))
end
abc_c = abc_smc(dist_c, [LogUniform(0.05, 5.0), LogUniform(20.0, 1000.0)]; n_particles = 100, generations = 5, rng = Xoshiro(6), names = [:h_max, :K], verbose = true)
save_csv("fig6c_particles.csv", ["h_max", "K", "weight"], hcat(abc_c.particles, abc_c.weights))
save_csv("fig6c_schedule.csv", ["generation", "epsilon", "acceptance"], hcat(1:length(abc_c.epsilons), abc_c.epsilons, abc_c.acceptance))
save_kv("fig6c_summary.csv", ["true_h_max" => truth_c.h_max, "true_K" => truth_c.K, "abc_mean_h_max" => posterior_mean(abc_c)[1], "abc_mean_K" => posterior_mean(abc_c)[2],
        "ci_h_lo" => credible_interval(abc_c, :h_max)[1], "ci_h_hi" => credible_interval(abc_c, :h_max)[2], "ci_K_lo" => credible_interval(abc_c, :K)[1], "ci_K_hi" => credible_interval(abc_c, :K)[2]])
# (d) posterior predictive kill curves from posterior-mean parameters versus the observed curve
pm = posterior_mean(abc_c)
pp = [simulate_population(rm_, xr, 400, (0.0, 100.0); settings = stc, perturbation = mkpert(pm[1], pm[2]), rng = Xoshiro(100 + s)) for s in 1:5]
i0 = argmin(abs.(obs_c.t .- 40.0))
save_csv("fig6d_ppc.csv", vcat(["t", "observed"], ["replicate_$s" for s in 1:5]),
         hcat(obs_c.t[i0:end], obs_c.popsize[i0:end] ./ obs_c.popsize[i0], [p.popsize[i0:end] ./ p.popsize[i0] for p in pp]...))
println("fig6 done")
