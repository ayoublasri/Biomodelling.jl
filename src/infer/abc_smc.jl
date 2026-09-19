# ---------------------------------------------------------------------------
# ABC-SMC (Toni et al. 2009) with adaptive tolerance schedule
# ---------------------------------------------------------------------------

"""
    ABCResult

Weighted particles (`particles`, n × d), `weights`, per-particle `distances`,
the tolerance schedule `epsilons`, per-generation `acceptance` rates and the
parameter `names`.
"""
struct ABCResult
    particles::Matrix{Float64}
    weights::Vector{Float64}
    distances::Vector{Float64}
    epsilons::Vector{Float64}
    acceptance::Vector{Float64}
    names::Vector{Symbol}
end
Base.show(io::IO, r::ABCResult) = print(io, "ABCResult(", size(r.particles, 1), " particles, ", size(r.particles, 2),
    " parameters, final ε=", round(r.epsilons[end], sigdigits = 3), ")")

"""
    abc_smc(distance, prior; n_particles=200, generations=6, alpha=0.5, rng, threads=true,
            names=Symbol[], min_acceptance=1e-3, verbose=false) -> ABCResult

Sequential Monte Carlo approximate Bayesian computation. `distance(θ, rng)` must
simulate with parameters `θ` (a vector ordered like `prior`) using the supplied
random number generator and return a distance to the observed summaries.
`prior` is a vector of univariate distributions. The tolerance of each
generation is the `alpha` quantile of the previous generation's distances.
"""
function abc_smc(distance::Function, prior::AbstractVector{<:UnivariateDistribution}; n_particles::Int=200,
                 generations::Int=6, alpha::Float64=0.5, rng::AbstractRNG=Random.default_rng(), threads::Bool=true,
                 names::Vector{Symbol}=Symbol[], min_acceptance::Float64=1e-3, verbose::Bool=false)
    d = length(prior)
    names = isempty(names) ? [Symbol("θ", i) for i in 1:d] : names
    logprior(θ) = sum(logpdf(prior[i], θ[i]) for i in 1:d)
    # generation 1: sample the prior
    P = zeros(n_particles, d)
    for i in 1:n_particles, k in 1:d
        P[i, k] = rand(rng, prior[k])
    end
    D = _eval_batch(distance, P, rng, threads)
    W = fill(1 / n_particles, n_particles)
    eps = [Inf]
    acc = [1.0]
    for gen in 2:generations
        ε = quantile(D, alpha)
        push!(eps, ε)
        σ = [sqrt(2 * max(var(view(P, :, k), Weights(W)), 1e-12)) for k in 1:d]
        newP = zeros(n_particles, d); newD = zeros(n_particles); newW = zeros(n_particles)
        nacc = 0; ntried = 0
        cum = cumsum(W)
        while nacc < n_particles
            batch = n_particles
            prop = zeros(batch, d)
            for b in 1:batch
                j = searchsortedfirst(cum, rand(rng) * cum[end])
                j = clamp(j, 1, n_particles)
                for k in 1:d
                    prop[b, k] = P[j, k] + σ[k] * randn(rng)
                end
            end
            valid = [isfinite(logprior(view(prop, b, :))) for b in 1:batch]
            dist = fill(Inf, batch)
            vidx = findall(valid)
            if !isempty(vidx)
                dist[vidx] = _eval_batch(distance, prop[vidx, :], rng, threads)
            end
            ntried += batch
            for b in 1:batch
                (valid[b] && dist[b] <= ε) || continue
                nacc == n_particles && break
                nacc += 1
                newP[nacc, :] = prop[b, :]
                newD[nacc] = dist[b]
                # importance weight
                den = 0.0
                for j in 1:n_particles
                    lk = 0.0
                    for k in 1:d
                        lk += logpdf(Normal(P[j, k], σ[k]), prop[b, k])
                    end
                    den += W[j] * exp(lk)
                end
                newW[nacc] = exp(logprior(view(prop, b, :))) / max(den, 1e-300)
            end
            if ntried > n_particles && nacc / ntried < min_acceptance
                @warn "ABC-SMC stopped at generation $gen: acceptance rate below $min_acceptance"
                return ABCResult(P, W, D, eps[1:end-1], acc, names)
            end
        end
        P = newP; D = newD; W = newW ./ sum(newW)
        push!(acc, n_particles / ntried)
        verbose && println("generation $gen: ε = $(round(ε, sigdigits=3)), acceptance = $(round(acc[end], sigdigits=3))")
    end
    ABCResult(P, W, D, eps, acc, names)
end

function _eval_batch(distance::Function, P::AbstractMatrix{Float64}, rng::AbstractRNG, threads::Bool)
    n = size(P, 1)
    seeds = rand(rng, UInt64, n)
    out = zeros(n)
    if threads && Threads.nthreads() > 1
        @sync for chunk in _chunks(n, Threads.nthreads())
            Threads.@spawn for i in chunk
                out[i] = Float64(distance(P[i, :], Xoshiro(seeds[i])))
            end
        end
    else
        for i in 1:n
            out[i] = Float64(distance(P[i, :], Xoshiro(seeds[i])))
        end
    end
    out
end

posterior_mean(r::ABCResult) = vec(sum(r.particles .* r.weights; dims = 1))

"""
    posterior_quantile(result, k, q) -> Float64
    posterior_median(result) -> Vector{Float64}

Weighted quantiles of the particle approximation (`k` is an index or a name).
"""
function posterior_quantile(r::ABCResult, k, q::Real)
    kk = k isa Integer ? Int(k) : findfirst(==(k), r.names)
    x = r.particles[:, kk]
    o = sortperm(x)
    cw = cumsum(r.weights[o]) ./ sum(r.weights)
    x[o[min(searchsortedfirst(cw, q), length(x))]]
end
posterior_median(r::ABCResult) = [posterior_quantile(r, k, 0.5) for k in 1:size(r.particles, 2)]

"""
    credible_interval(result, k; level=0.95) -> (lo, hi)

Weighted quantile interval of parameter `k` (index or name).
"""
function credible_interval(r::ABCResult, k; level::Float64=0.95)
    kk = k isa Integer ? Int(k) : findfirst(==(k), r.names)
    x = r.particles[:, kk]
    o = sortperm(x)
    cw = cumsum(r.weights[o]) ./ sum(r.weights)
    lo = x[o[searchsortedfirst(cw, (1 - level) / 2)]]
    hi = x[o[min(searchsortedfirst(cw, 1 - (1 - level) / 2), length(x))]]
    (lo, hi)
end
