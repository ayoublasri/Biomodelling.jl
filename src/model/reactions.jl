# ---------------------------------------------------------------------------
# Reaction and kinetics specification (user-facing, symbolic)
# ---------------------------------------------------------------------------

"""
    Kinetics

Abstract supertype of propensity laws. Concrete types: [`MassAction`](@ref),
[`Hill`](@ref), [`Custom`](@ref).
"""
abstract type Kinetics end

"""
    MassAction(rate::Symbol)

Mass-action kinetics. The propensity is `k * prod_i binomial(x_i, ν_i) * V^e`,
where `k` is the named parameter `rate`, `x_i` and `ν_i` are the reactant counts
and stoichiometric coefficients, and `e` is the volume exponent of the reaction
(see [`Reaction`](@ref)).
"""
struct MassAction <: Kinetics
    rate::Symbol
end

const ParamOrValue = Union{Symbol,Float64}
_pv(x::Symbol) = x
_pv(x::Real) = Float64(x)
function _expand(v, m::Int)
    if v isa AbstractVector
        length(v) == m || throw(ArgumentError("expected $m values, got $(length(v))"))
        return ParamOrValue[_pv(x) for x in v]
    end
    return ParamOrValue[_pv(v) for _ in 1:m]
end

"""
    Hill(rate::Symbol; activators=[], inhibitors=[], K, n=2.0, basal=0.0,
         logic=:and, concentration=true)

Regulated kinetics. The propensity is

    k * V^e * [basal + (1 - basal) * F(c)] * prod_i binomial(x_i, ν_i)

where `F` combines one Hill term per regulator, `h_act(c) = c^n / (K^n + c^n)` for
activators and `h_inh(c) = K^n / (K^n + c^n)` for inhibitors, with `c` the
regulator concentration (`x / V`) when `concentration = true`, or its count
otherwise. `logic = :and` multiplies the terms (all regulators needed),
`logic = :or` combines them as `1 - prod(1 - h_i)`.

`K`, `n` and `basal` may be parameter names (`Symbol`) or numeric literals; a
scalar is broadcast to all regulators, a vector gives one value per regulator
(activators first, then inhibitors). Literals become anonymous model parameters.
"""
struct Hill <: Kinetics
    rate::Symbol
    regulators::Vector{Symbol}
    modes::Vector{Symbol}
    K::Vector{ParamOrValue}
    n::Vector{ParamOrValue}
    basal::ParamOrValue
    logic::Symbol
    concentration::Bool
end

function Hill(rate::Symbol; activators=Symbol[], inhibitors=Symbol[], K, n=2.0,
              basal=0.0, logic::Symbol=:and, concentration::Bool=true)
    acts = activators isa Symbol ? Symbol[activators] : collect(Symbol, activators)
    inhs = inhibitors isa Symbol ? Symbol[inhibitors] : collect(Symbol, inhibitors)
    regs = vcat(acts, inhs)
    isempty(regs) && throw(ArgumentError("Hill kinetics need at least one regulator"))
    logic in (:and, :or) || throw(ArgumentError("logic must be :and or :or"))
    modes = vcat(fill(:activate, length(acts)), fill(:inhibit, length(inhs)))
    Hill(rate, regs, modes, _expand(K, length(regs)), _expand(n, length(regs)),
         _pv(basal), logic, concentration)
end

"""
    Custom(f; params=Symbol[])

User-defined propensity `f(x, V, θ, t)`, where `x` is the vector of counts, `V`
the cell volume, `θ` the values of the named parameters listed in `params` (in
that order) and `t` the time. The propensity is multiplied by the mass-action
reactant factor and the reaction's volume exponent, like the other kinetics.
Custom propensities are assumed to depend on every species.
"""
struct Custom <: Kinetics
    f::Function
    params::Vector{Symbol}
end
Custom(f; params=Symbol[]) = Custom(f, collect(Symbol, params))

"""
    Reaction(name, reactants, products, kinetics; volume=:auto, copy_number=false)

A reaction channel. `reactants` and `products` are vectors of species symbols
(coefficient 1) or `Symbol => coefficient` pairs; `:NULL` and `:∅` are ignored.

`volume` sets how the propensity scales with cell volume `V`:
`:auto` uses `V^(1 - order)` (zero-order production ∝ V, first order unchanged,
bimolecular ∝ 1/V), `:none`, `:proportional` (∝ V) and `:inverse` (∝ 1/V).

`copy_number = true` multiplies the propensity by the cell's gene copy number
(1 before replication, 2 after); use it for constitutive transcription of genes
without an explicit promoter.
"""
struct Reaction
    name::String
    reactants::Vector{Pair{Symbol,Int}}
    products::Vector{Pair{Symbol,Int}}
    kinetics::Kinetics
    volume::Symbol
    copy_number::Bool
end

const _NULLS = (:NULL, :∅, :nothing, :none)
_topairs(::Nothing) = Pair{Symbol,Int}[]
_topairs(v::Pair) = Pair{Symbol,Int}[Symbol(first(v)) => Int(last(v))]
_topairs(s::Symbol) = s in _NULLS ? Pair{Symbol,Int}[] : Pair{Symbol,Int}[s => 1]
_topairs(t::Tuple) = _topairs(collect(t))
function _topairs(v::AbstractVector)
    out = Pair{Symbol,Int}[]
    for x in v
        if x isa Pair
            push!(out, Symbol(first(x)) => Int(last(x)))
        elseif x isa Symbol
            x in _NULLS || push!(out, x => 1)
        else
            throw(ArgumentError("cannot interpret $x as a species"))
        end
    end
    out
end

function Reaction(name::AbstractString, reactants, products, kinetics::Kinetics;
                  volume::Symbol=:auto, copy_number::Bool=false)
    volume in (:auto, :none, :proportional, :inverse) ||
        throw(ArgumentError("volume must be :auto, :none, :proportional or :inverse"))
    Reaction(String(name), _topairs(reactants), _topairs(products), kinetics, volume, copy_number)
end

reaction_order(r::Reaction) = sum(last, r.reactants; init=0)
function Base.show(io::IO, r::Reaction)
    lhs = isempty(r.reactants) ? "∅" : join(("$(c == 1 ? "" : "$c ")$s" for (s, c) in r.reactants), " + ")
    rhs = isempty(r.products) ? "∅" : join(("$(c == 1 ? "" : "$c ")$s" for (s, c) in r.products), " + ")
    print(io, "Reaction(\"", r.name, "\": ", lhs, " → ", rhs, ", ", nameof(typeof(r.kinetics)), ")")
end
