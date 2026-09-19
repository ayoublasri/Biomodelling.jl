# ---------------------------------------------------------------------------
# Exact stationary distribution of the telegraph model (Peccoud and Ycart 1995)
# ---------------------------------------------------------------------------

"""
    kummer_1f1_neg(a, b, λ) -> log 1F1(a; b; -λ)

Log of Kummer's confluent hypergeometric function at negative argument, via
Kummer's transformation `1F1(a;b;-λ) = e^{-λ} 1F1(b-a; b; λ)`, whose series has
positive terms. Partial sums are rescaled to avoid overflow for large `λ`.
"""
function kummer_1f1_neg(a::Float64, b::Float64, λ::Float64; maxterms::Int=2_000_000, rtol::Float64=1e-14)
    ap = b - a
    term = 1.0
    s = 1.0
    acc = 0.0
    k = 0
    while k < maxterms
        term *= (ap + k) / (b + k) * λ / (k + 1)
        s += term
        k += 1
        if s > 1e200
            s *= 1e-200; term *= 1e-200; acc += 200 * log(10.0)
        end
        (k > λ && term < rtol * s) && break
    end
    log(s) + acc - λ
end

"""
    telegraph_pmf(n, k_on, k_off, k_tx, γ=1.0) -> Float64
    telegraph_pmf(0:nmax, ...) -> Vector{Float64}

Stationary probability of `n` mRNA molecules in the two-state telegraph model
with activation `k_on`, inactivation `k_off`, transcription `k_tx` (from the
active state) and degradation `γ` (Beta-Poisson distribution, Peccoud and Ycart
1995).
"""
function telegraph_pmf(n::Integer, k_on::Real, k_off::Real, k_tx::Real, γ::Real=1.0)
    a = k_on / γ; b = k_off / γ; λ = k_tx / γ
    (a > 0 && b > 0 && λ >= 0) || throw(ArgumentError("rates must be positive"))
    λ == 0 && return n == 0 ? 1.0 : 0.0
    logp = n * log(λ) - loggamma(n + 1) + loggamma(a + n) + loggamma(a + b) - loggamma(a + b + n) - loggamma(a) +
           kummer_1f1_neg(a + n, a + b + n, λ)
    p = exp(logp)
    isfinite(p) ? min(p, 1.0) : 0.0
end
telegraph_pmf(ns::AbstractVector{<:Integer}, k_on, k_off, k_tx, γ=1.0) = [telegraph_pmf(n, k_on, k_off, k_tx, γ) for n in ns]

"""
    telegraph_loglik(counts, k_on, k_off, k_tx; γ=1.0) -> Float64
"""
function telegraph_loglik(counts::AbstractVector{<:Integer}, k_on::Real, k_off::Real, k_tx::Real; γ::Real=1.0)
    cm = countmap(counts)
    ll = 0.0
    for (n, c) in cm
        p = telegraph_pmf(n, k_on, k_off, k_tx, γ)
        ll += c * log(max(p, 1e-300))
    end
    ll
end

"""
    fit_telegraph(counts; γ=1.0, init=nothing) -> NamedTuple

Maximum-likelihood telegraph parameters from a sample of mRNA counts (rates
returned in units of `γ`). A two-dimensional grid search with the transcription
rate tied to the sample mean initialises two rounds of Nelder-Mead optimisation
in log-parameter space.
"""
function fit_telegraph(counts::AbstractVector{<:Integer}; γ::Real=1.0, init=nothing)
    m = mean(counts)
    m <= 0 && throw(ArgumentError("counts must have a positive mean"))
    λmax = 20.0 * (maximum(counts) + 10)
    function nll(θ)
        a, b, λ = exp(θ[1]), exp(θ[2]), exp(θ[3])
        (λ > λmax || a > 1e4 || b > 1e4 || a < 1e-6 || b < 1e-6) && return Inf
        v = -telegraph_loglik(counts, a, b, λ)
        isfinite(v) ? v : Inf
    end
    θ0 = if init === nothing
        best = (Inf, zeros(3))
        for la in range(-5, 4; length = 19), lb in range(-5, 4; length = 19)
            a = exp(la); b = exp(lb)
            θ = [la, lb, log(m * (a + b) / a)]
            v = nll(θ)
            v < best[1] && (best = (v, θ))
        end
        best[2]
    else
        log.(collect(Float64, init))
    end
    r1 = Optim.optimize(nll, θ0, Optim.NelderMead(), Optim.Options(iterations = 4000))
    r2 = Optim.optimize(nll, Optim.minimizer(r1), Optim.NelderMead(), Optim.Options(iterations = 4000, g_tol = 1e-10))
    θ = Optim.minimizer(r2)
    (k_on = exp(θ[1]) * γ, k_off = exp(θ[2]) * γ, k_tx = exp(θ[3]) * γ, loglik = -Optim.minimum(r2),
     converged = Optim.converged(r2) || Optim.converged(r1))
end
