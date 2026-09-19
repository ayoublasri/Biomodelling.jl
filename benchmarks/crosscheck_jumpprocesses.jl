# Cross-check of Biomodelling's DirectSSA against JumpProcesses.jl (SciML) on the
# telegraph model. Run in an environment containing JumpProcesses.
using Biomodelling, JumpProcesses, Random, Statistics, Distributions
using Random: Xoshiro

k_on, k_off, k_tx, k_dm = 0.4, 0.6, 12.0, 1.0
T, n = 30.0, 20_000

# --- Biomodelling ---------------------------------------------------------------
tm = telegraph_model(; k_on, k_off, k_tx, k_dm, volume_scaled = false)
x0 = initial_state(tm; G_off = 1)
tb = @elapsed Eb = ensemble_final(tm, x0, T, n; rng = Xoshiro(1))
cb = Eb[:, speciesindex(tm, :mRNA)]

# --- JumpProcesses ------------------------------------------------------------
# species order: G_off, G_on, mRNA
rate1(u, p, t) = p[1] * u[1]; affect1!(i) = (i.u[1] -= 1; i.u[2] += 1; nothing)
rate2(u, p, t) = p[2] * u[2]; affect2!(i) = (i.u[2] -= 1; i.u[1] += 1; nothing)
rate3(u, p, t) = p[3] * u[2]; affect3!(i) = (i.u[3] += 1; nothing)
rate4(u, p, t) = p[4] * u[3]; affect4!(i) = (i.u[3] -= 1; nothing)
jumps = (ConstantRateJump(rate1, affect1!), ConstantRateJump(rate2, affect2!),
         ConstantRateJump(rate3, affect3!), ConstantRateJump(rate4, affect4!))
dprob = DiscreteProblem([1, 0, 0], (0.0, T), [k_on, k_off, k_tx, k_dm])
jprob = JumpProblem(dprob, Direct(), jumps...; save_positions = (false, false), rng = Xoshiro(2))
eprob = EnsembleProblem(jprob)
tj = @elapsed sol = solve(eprob, SSAStepper(), EnsembleSerial(); trajectories = n, saveat = [T])
cj = [Int(s.u[end][3]) for s in sol]

# --- comparison ----------------------------------------------------------------
pmf = telegraph_pmf(0:100, k_on, k_off, k_tx, k_dm)
ks(sample) = maximum(abs.([count(<=(k), sample) / length(sample) for k in 0:100] .- cumsum(pmf)))
println("Biomodelling : mean=$(round(mean(cb), digits=3)) var=$(round(var(cb), digits=3)) KS=$(round(ks(cb), sigdigits=3)) time=$(round(tb, digits=2))s")
println("JumpProcesses: mean=$(round(mean(cj), digits=3)) var=$(round(var(cj), digits=3)) KS=$(round(ks(cj), sigdigits=3)) time=$(round(tj, digits=2))s")
println("theory       : mean=$(round(k_tx*k_on/(k_on+k_off)/k_dm, digits=3)) var=$(round(sum(pmf .* (0:100).^2) - sum(pmf .* (0:100))^2, digits=3))")
# two-sample KS between the two simulators
emp(sample) = [count(<=(k), sample) / length(sample) for k in 0:100]
println("two-sample KS(Biomodelling, JumpProcesses) = ", round(maximum(abs.(emp(cb) .- emp(cj))), sigdigits = 3),
        "  (99% critical ≈ ", round(1.63 * sqrt(2 / n), sigdigits = 3), ")")
