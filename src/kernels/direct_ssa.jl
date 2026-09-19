"""
    DirectSSA()

Gillespie's direct method (exact). Propensities are updated through the model's
dependency graph after every event.
"""
struct DirectSSA <: AbstractKernel end

"""
    simulate!(x, model, p, V, copies, t0, t1, kernel, ws, rng) -> t1

Advance the state `x` in place from `t0` to `t1` at fixed volume `V` and gene
copy number `copies`.
"""
function simulate!(x::AbstractVector{Int}, m::ReactionModel, p::AbstractVector{Float64}, V::Float64,
                   copies::Int, t0::Float64, t1::Float64, ::DirectSSA, ws::Workspace, rng::AbstractRNG)
    a = ws.a
    a0 = propensities!(a, m, x, V, copies, p, t0)
    t = t0
    nev = 0
    M = length(a)
    while true
        a0 <= 0.0 && return t1
        t += randexp(rng) / a0
        t >= t1 && return t1
        r = rand(rng) * a0
        j = 1
        s = a[1]
        @inbounds while s < r && j < M
            j += 1
            s += a[j]
        end
        fire!(x, m, j)
        a0 = update_propensities!(a, m, j, x, V, copies, p, t, a0)
        nev += 1
        if nev % 512 == 0
            a0 = sum(a)
        elseif a0 < 0.0
            a0 = sum(a)
        end
    end
end

"""
    ssa_events!(x, model, p, V, copies, t0, t1, nevents, ws, rng) -> t

Perform at most `nevents` exact events, stopping at `t1`. Returns the time reached.
"""
function ssa_events!(x::AbstractVector{Int}, m::ReactionModel, p::AbstractVector{Float64}, V::Float64,
                     copies::Int, t0::Float64, t1::Float64, nevents::Int, ws::Workspace, rng::AbstractRNG)
    a = ws.a
    a0 = propensities!(a, m, x, V, copies, p, t0)
    t = t0
    M = length(a)
    for _ in 1:nevents
        a0 <= 0.0 && return t1
        τ = randexp(rng) / a0
        t + τ >= t1 && return t1
        t += τ
        r = rand(rng) * a0
        j = 1
        s = a[1]
        @inbounds while s < r && j < M
            j += 1
            s += a[j]
        end
        fire!(x, m, j)
        a0 = update_propensities!(a, m, j, x, V, copies, p, t, a0)
        a0 < 0.0 && (a0 = sum(a))
    end
    t
end
