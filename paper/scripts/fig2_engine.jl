# Figure 2: engine validation against exact results, between kernels, against JumpProcesses.jl, and runtime.
include(joinpath(@__DIR__, "common.jl"))
using JumpProcesses

# (a) birth-death → Poisson
bd = birth_death_model(k = 12.0, γ = 1.0)
E = ensemble_final(bd, [12], 15.0, 50_000; rng = Xoshiro(1))[:, 1]
save_csv("fig2a_birthdeath.csv", ["n", "empirical", "poisson"], hcat(0:40, empirical_pmf(E, 40), pdf.(Poisson(12.0), 0:40)))

# (b) telegraph → Beta-Poisson (two regimes)
rows = Any[]
for (k_on, k_off, k_tx, tag) in ((0.4, 0.6, 12.0, "intermediate"), (0.05, 0.5, 40.0, "bursty"))
    tm = telegraph_model(; k_on, k_off, k_tx, k_dm = 1.0, volume_scaled = false)
    E = ensemble_final(tm, initial_state(tm; G_off = 1), 40.0, 50_000; rng = Xoshiro(2))[:, speciesindex(tm, :mRNA)]
    pm = telegraph_pmf(0:80, k_on, k_off, k_tx, 1.0)
    emp = empirical_pmf(E, 80)
    for n in 0:80
        push!(rows, [tag, n, emp[n+1], pm[n+1]])
    end
end
save_csv("fig2b_telegraph.csv", ["regime", "n", "empirical", "beta_poisson"], permutedims(reduce(hcat, rows)))

# (c) bursty protein → negative binomial
a, b, γ = 1.5, 6.0, 1.0
bp = bursty_protein_model(; a, b, γ)
E = ensemble_final(bp, initial_state(bp), 25.0, 50_000; rng = Xoshiro(3))[:, speciesindex(bp, :protein)]
nb = NegativeBinomial(a / γ, 1 / (1 + b))
save_csv("fig2c_bursty.csv", ["n", "empirical", "negbin"], hcat(0:80, empirical_pmf(E, 80), pdf.(nb, 0:80)))

# (d) kernel agreement: KS distance to the exact law, plus runtime per 10^4 cells
rows = Any[]
tm = telegraph_model(k_on = 0.4, k_off = 0.6, k_tx = 12.0, k_dm = 1.0; volume_scaled = false)
x0 = initial_state(tm; G_off = 1)
pm = telegraph_pmf(0:100, 0.4, 0.6, 12.0, 1.0)
pmbd = pdf.(Poisson(12.0), 0:100)
for (name, kern) in (("DirectSSA", DirectSSA()), ("HybridSSATau(0.05)", HybridSSATau(0.05)), ("HybridSSATau(0.2)", HybridSSATau(0.2)),
                     ("AdaptiveTauLeap(0.03)", AdaptiveTauLeap(ε = 0.03)), ("AdaptiveTauLeap(0.1)", AdaptiveTauLeap(ε = 0.1)))
    ensemble_final(tm, x0, 1.0, 100; kernel = kern, rng = Xoshiro(4))   # compile
    t = @elapsed E = ensemble_final(tm, x0, 40.0, 20_000; kernel = kern, rng = Xoshiro(4))
    t2 = @elapsed E2 = ensemble_final(bd, [12], 15.0, 20_000; kernel = kern, rng = Xoshiro(5))
    push!(rows, [name, ks_discrete(E[:, 3], pm), t / 2, ks_discrete(E2[:, 1], pmbd), t2 / 2])
end
save_csv("fig2d_kernels.csv", ["kernel", "ks_telegraph", "sec_per_1e4_telegraph", "ks_birthdeath", "sec_per_1e4_birthdeath"], permutedims(reduce(hcat, rows)))

# (e) cross-check against JumpProcesses.jl (SciML) on the telegraph model
n = 20_000
rate1(u, p, t) = p[1] * u[1]; affect1!(i) = (i.u[1] -= 1; i.u[2] += 1; nothing)
rate2(u, p, t) = p[2] * u[2]; affect2!(i) = (i.u[2] -= 1; i.u[1] += 1; nothing)
rate3(u, p, t) = p[3] * u[2]; affect3!(i) = (i.u[3] += 1; nothing)
rate4(u, p, t) = p[4] * u[3]; affect4!(i) = (i.u[3] -= 1; nothing)
jumps = (ConstantRateJump(rate1, affect1!), ConstantRateJump(rate2, affect2!), ConstantRateJump(rate3, affect3!), ConstantRateJump(rate4, affect4!))
dprob = DiscreteProblem([1, 0, 0], (0.0, 40.0), [0.4, 0.6, 12.0, 1.0])
jprob = JumpProblem(dprob, Direct(), jumps...; save_positions = (false, false), rng = Xoshiro(6))
solve(EnsembleProblem(jprob), SSAStepper(), EnsembleSerial(); trajectories = 10, saveat = [40.0])
tj = @elapsed sol = solve(EnsembleProblem(jprob), SSAStepper(), EnsembleSerial(); trajectories = n, saveat = [40.0])
cj = [Int(sol.u[i].u[end][3]) for i in eachindex(sol.u)]   # EnsembleSolution iterates over states, not trajectories
tb = @elapsed cb = ensemble_final(tm, x0, 40.0, n; rng = Xoshiro(7), threads = false)[:, 3]
save_csv("fig2e_crosscheck.csv", ["n", "biomodelling", "jumpprocesses", "exact"], hcat(0:80, empirical_pmf(cb, 80), empirical_pmf(cj, 80), pm[1:81]))
save_kv("fig2e_crosscheck_stats.csv", ["ks_biomodelling_exact" => ks_discrete(cb, pm), "ks_jumpprocesses_exact" => ks_discrete(cj, pm),
        "ks_two_sample" => maximum(abs.(empirical_pmf(cb, 100) .- empirical_pmf(cj, 100))), "ks_two_sample_crit99" => 1.63 * sqrt(2 / n),
        "time_biomodelling_serial_s" => tb, "time_jumpprocesses_serial_s" => tj, "mean_biomodelling" => mean(cb), "mean_jumpprocesses" => mean(cj),
        "mean_exact" => sum(pm .* (0:100))])

# (f) runtime scaling of population simulations
rows = Any[]
for (ng, nc) in ((10, 100), (10, 1000), (10, 5000), (50, 1000), (100, 1000)), (kname, k) in (("DirectSSA", DirectSSA()), ("HybridSSATau", HybridSSATau(0.1)))
    m, _ = random_grn(ng; n_activations = ng, n_inhibitions = ng ÷ 2, telegraph = true, rng = Xoshiro(8))
    xg = zeros(Int, nspecies(m)); for g in m.promoter_groups; xg[g[1]] = 1; end
    st = PopulationSettings(dt = 0.1, kernel = k, growth = ExponentialGrowth(log(2) / 20), track_lineage = false, record_every = 20)
    simulate_population(m, xg, nc, (0.0, 1.0); settings = st, rng = Xoshiro(9))
    t = @elapsed simulate_population(m, xg, nc, (0.0, 20.0); settings = st, rng = Xoshiro(9))
    push!(rows, [ng, nc, kname, Threads.nthreads(), t, nc * 200 / t])
    @printf("%4d genes %5d cells %-14s %7.2fs\n", ng, nc, kname, t)
end
save_csv("fig2f_runtime.csv", ["genes", "cells", "kernel", "threads", "seconds", "cell_steps_per_second"], permutedims(reduce(hcat, rows)))
println("fig2 done")
