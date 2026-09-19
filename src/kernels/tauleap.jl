"""
    TauLeap(τ)

Fixed-step Poisson tau-leaping. Throws an error if a leap would drive a species
negative; use [`HybridSSATau`](@ref) for automatic fallback.
"""
struct TauLeap <: AbstractKernel
    τ::Float64
end

"""
    HybridSSATau(τ)

Fixed-step Poisson tau-leaping that falls back to the exact direct method for any
leap interval in which a species would become negative (the "switching" scheme
of the original package, made exact over the rejected interval).
"""
struct HybridSSATau <: AbstractKernel
    τ::Float64
end

function _leap!(xtmp::Vector{Int}, x::AbstractVector{Int}, m::ReactionModel, a::Vector{Float64},
                τ::Float64, rng::AbstractRNG)
    copyto!(xtmp, x)
    @inbounds for j in eachindex(a)
        aj = a[j]
        aj <= 0.0 && continue
        n = pois_rand(rng, aj * τ)
        n == 0 && continue
        for (s, c) in m.stoich[j]
            xtmp[s] += n * c
        end
    end
    @inbounds for s in eachindex(xtmp)
        xtmp[s] < 0 && return false
    end
    true
end

function simulate!(x::AbstractVector{Int}, m::ReactionModel, p::AbstractVector{Float64}, V::Float64,
                   copies::Int, t0::Float64, t1::Float64, k::Union{TauLeap,HybridSSATau}, ws::Workspace, rng::AbstractRNG)
    t = t0
    a = ws.a
    while t < t1 - 1e-12 * max(1.0, abs(t1))
        τ = min(k.τ, t1 - t)
        propensities!(a, m, x, V, copies, p, t)
        ok = _leap!(ws.xtmp, x, m, a, τ, rng)
        if ok
            copyto!(x, ws.xtmp)
        elseif k isa TauLeap
            error("tau-leap produced a negative population at t=$t; use a smaller τ or HybridSSATau")
        else
            simulate!(x, m, p, V, copies, t, t + τ, DirectSSA(), ws, rng)
        end
        t += τ
    end
    t1
end
