# ---------------------------------------------------------------------------
# Propensity evaluation
# ---------------------------------------------------------------------------

@inline function _combos(x::Int, ν::Int)
    ν == 1 && return float(x)
    x < ν && return 0.0
    ν == 2 && return x * (x - 1) / 2
    ν == 3 && return x * (x - 1) * (x - 2) / 6
    return float(binomial(x, ν))
end

@inline function _massfactor(x::AbstractVector{Int}, reactants::Vector{Tuple{Int,Int}})
    f = 1.0
    @inbounds for (s, ν) in reactants
        f *= _combos(x[s], ν)
        f == 0.0 && return 0.0
    end
    f
end

@inline function _volfactor(V::Float64, e::Float64)
    e == 0.0 && return 1.0
    e == 1.0 && return V
    e == -1.0 && return 1.0 / V
    return V^e
end

@inline function _hill(c::Float64, K::Float64, n::Float64)
    c <= 0.0 && return 0.0
    cn = c^n
    return cn / (K^n + cn)
end

@inline function _regfactor(k::HillC, x::AbstractVector{Int}, V::Float64, p::AbstractVector{Float64})
    f = 1.0
    if k.and_logic
        @inbounds for i in eachindex(k.regs)
            c = k.conc ? x[k.regs[i]] / V : float(x[k.regs[i]])
            h = _hill(c, p[k.K[i]], p[k.n[i]])
            f *= k.activate[i] ? h : 1.0 - h
        end
    else
        g = 1.0
        @inbounds for i in eachindex(k.regs)
            c = k.conc ? x[k.regs[i]] / V : float(x[k.regs[i]])
            h = _hill(c, p[k.K[i]], p[k.n[i]])
            g *= 1.0 - (k.activate[i] ? h : 1.0 - h)
        end
        f = 1.0 - g
    end
    b = p[k.basal]
    return b + (1.0 - b) * f
end

@inline function propensity(k::MassActionC, x::AbstractVector{Int}, V::Float64, copies::Int,
                            p::AbstractVector{Float64}, t::Float64)
    a = p[k.rate] * _massfactor(x, k.reactants) * _volfactor(V, k.vol_exp)
    k.copy_number ? a * copies : a
end

@inline function propensity(k::HillC, x::AbstractVector{Int}, V::Float64, copies::Int,
                            p::AbstractVector{Float64}, t::Float64)
    a = p[k.rate] * _massfactor(x, k.reactants) * _volfactor(V, k.vol_exp) * _regfactor(k, x, V, p)
    k.copy_number ? a * copies : a
end

function propensity(k::CustomC, x::AbstractVector{Int}, V::Float64, copies::Int,
                    p::AbstractVector{Float64}, t::Float64)
    θ = view(p, k.params)
    a = float(k.f(x, V, θ, t)) * _massfactor(x, k.reactants) * _volfactor(V, k.vol_exp)
    k.copy_number ? a * copies : a
end

"""
    propensities!(a, model, x, V, copies, p, t) -> a0

Evaluate all propensities into `a` and return their sum.
"""
function propensities!(a::Vector{Float64}, m::ReactionModel, x::AbstractVector{Int}, V::Float64,
                       copies::Int, p::AbstractVector{Float64}, t::Float64)
    a0 = 0.0
    @inbounds for j in eachindex(a)
        aj = propensity(m.kinetics[j], x, V, copies, p, t)
        a[j] = aj
        a0 += aj
    end
    a0
end

"""
    propensities(model, x; V=1.0, copies=1, p=model.p0, t=0.0) -> Vector{Float64}
"""
function propensities(m::ReactionModel, x::AbstractVector{<:Integer}; V::Real=1.0, copies::Integer=1,
                      p::AbstractVector{<:Real}=m.p0, t::Real=0.0)
    a = zeros(nreactions(m))
    propensities!(a, m, Vector{Int}(x), Float64(V), Int(copies), Vector{Float64}(p), Float64(t))
    a
end

@inline function update_propensities!(a::Vector{Float64}, m::ReactionModel, j::Int, x::AbstractVector{Int},
                                      V::Float64, copies::Int, p::AbstractVector{Float64}, t::Float64, a0::Float64)
    @inbounds for k in m.depgraph[j]
        old = a[k]
        new = propensity(m.kinetics[k], x, V, copies, p, t)
        a[k] = new
        a0 += new - old
    end
    a0
end

@inline function fire!(x::AbstractVector{Int}, m::ReactionModel, j::Int, n::Int=1)
    @inbounds for (s, c) in m.stoich[j]
        x[s] += n * c
    end
    x
end
