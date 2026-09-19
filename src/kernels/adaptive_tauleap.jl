"""
    AdaptiveTauLeap(; ε=0.03, n_critical=10, ssa_multiple=10.0, ssa_events=100)

Adaptive tau-leaping with the step-size selection of Cao, Gillespie and Petzold
(2006) and the critical-reaction partition of Cao, Gillespie and Petzold (2005).
Reactions within `n_critical` firings of exhausting a reactant are simulated
exactly; when the proposed leap is shorter than `ssa_multiple / a0`, `ssa_events`
exact events are performed instead.
"""
Base.@kwdef struct AdaptiveTauLeap <: AbstractKernel
    ε::Float64 = 0.03
    n_critical::Int = 10
    ssa_multiple::Float64 = 10.0
    ssa_events::Int = 100
end

@inline function _gfactor(hor::Int, mrr::Int, x::Int)
    hor <= 1 && return 1.0
    if hor == 2
        mrr == 1 && return 2.0
        return x > 1 ? 2.0 + 1.0 / (x - 1) : 2.0
    end
    # hor >= 3
    mrr == 1 && return 3.0
    mrr == 2 && return x > 1 ? 1.5 * (2.0 + 1.0 / (x - 1)) : 3.0
    return x > 2 ? 3.0 + 1.0 / (x - 1) + 2.0 / (x - 2) : 3.0
end

function _tau_noncritical(m::ReactionModel, x::AbstractVector{Int}, a::Vector{Float64}, crit::Vector{Bool},
                          ε::Float64, mu::Vector{Float64}, sig::Vector{Float64})
    fill!(mu, 0.0); fill!(sig, 0.0)
    @inbounds for j in eachindex(a)
        (crit[j] || a[j] <= 0.0) && continue
        for (s, c) in m.stoich[j]
            mu[s] += c * a[j]
            sig[s] += c * c * a[j]
        end
    end
    τ = Inf
    @inbounds for s in eachindex(mu)
        m.hor[s] == 0 && continue
        h = max(ε * x[s] / _gfactor(m.hor[s], m.mrr[s], x[s]), 1.0)
        mu[s] != 0.0 && (τ = min(τ, h / abs(mu[s])))
        sig[s] != 0.0 && (τ = min(τ, h * h / sig[s]))
    end
    τ
end

function simulate!(x::AbstractVector{Int}, m::ReactionModel, p::AbstractVector{Float64}, V::Float64,
                   copies::Int, t0::Float64, t1::Float64, k::AdaptiveTauLeap, ws::Workspace, rng::AbstractRNG)
    t = t0
    a = ws.a
    crit = ws.critical
    while t < t1 - 1e-12 * max(1.0, abs(t1))
        a0 = propensities!(a, m, x, V, copies, p, t)
        a0 <= 0.0 && return t1
        # critical reactions
        a0c = 0.0
        @inbounds for j in eachindex(a)
            c = false
            if a[j] > 0.0
                L = typemax(Int)
                for (s, ν) in m.reactant_stoich[j]
                    L = min(L, fld(x[s], ν))
                end
                c = L < k.n_critical
            end
            crit[j] = c
            c && (a0c += a[j])
        end
        τ1 = _tau_noncritical(m, x, a, crit, k.ε, ws.mu, ws.sig)
        if τ1 < k.ssa_multiple / a0
            t = ssa_events!(x, m, p, V, copies, t, t1, k.ssa_events, ws, rng)
            continue
        end
        τ2 = a0c > 0.0 ? randexp(rng) / a0c : Inf
        τ = min(τ1, τ2, t1 - t)
        fire_critical = (τ2 <= τ1) && (τ2 <= t1 - t)
        accepted = false
        while !accepted
            copyto!(ws.xtmp, x)
            @inbounds for j in eachindex(a)
                (crit[j] || a[j] <= 0.0) && continue
                n = pois_rand(rng, a[j] * τ)
                n == 0 && continue
                for (s, c) in m.stoich[j]
                    ws.xtmp[s] += n * c
                end
            end
            if fire_critical
                r = rand(rng) * a0c
                s = 0.0
                jc = 0
                @inbounds for j in eachindex(a)
                    crit[j] || continue
                    s += a[j]
                    jc = j
                    s >= r && break
                end
                jc > 0 && fire!(ws.xtmp, m, jc)
            end
            accepted = all(>=(0), ws.xtmp)
            if !accepted
                τ /= 2
                fire_critical = false
            end
        end
        copyto!(x, ws.xtmp)
        t += τ
    end
    t1
end
