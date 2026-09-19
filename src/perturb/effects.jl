abstract type DrugEffect end

"""
    DeathHazard(; h_max, EC50, m=2.0, protect=nothing, K=1.0, q=2.0, concentration=true)

Dose-dependent death hazard `h = h_max d^m / (EC50^m + d^m)`, optionally reduced
by a protective species (`protect`, e.g. a resistance protein) through the
factor `K^q / (K^q + c^q)` where `c` is that species' concentration (or count if
`concentration = false`). This is the state-dependent survival mechanism of
Lasri and Sturrock (2020).
"""
struct DeathHazard <: DrugEffect
    h_max::Float64
    EC50::Float64
    m::Float64
    protect::Union{Nothing,Symbol}
    K::Float64
    q::Float64
    concentration::Bool
end
DeathHazard(; h_max, EC50, m=2.0, protect=nothing, K=1.0, q=2.0, concentration=true) =
    DeathHazard(Float64(h_max), Float64(EC50), Float64(m), protect, Float64(K), Float64(q), concentration)

"""
    GrowthInhibition(; IC50, m=2.0)

Multiplies the growth rate by `1 / (1 + (d / IC50)^m)`.
"""
struct GrowthInhibition <: DrugEffect
    IC50::Float64
    m::Float64
end
GrowthInhibition(; IC50, m=2.0) = GrowthInhibition(Float64(IC50), Float64(m))

"""
    RateModulation(param, f)

Multiplies the model parameter `param` by `f(d)` at dose `d` (e.g. drug-induced
transcription, or a memory modulator that increases promoter switching rates).
"""
struct RateModulation <: DrugEffect
    param::Symbol
    f::Function
end

"""
    GenePerturbation(param, factor; fraction=1.0, t_start=0.0, t_end=Inf)

Genetic perturbation (knockdown `factor < 1`, overexpression `factor > 1`,
knockout `factor = 0`) of parameter `param` in a random `fraction` of the cells
from `t_start` to `t_end`. The perturbed flag is inherited by daughters.
"""
struct GenePerturbation
    param::Symbol
    factor::Float64
    fraction::Float64
    t_start::Float64
    t_end::Float64
end
GenePerturbation(param::Symbol, factor::Real; fraction::Real=1.0, t_start::Real=0.0, t_end::Real=Inf) =
    GenePerturbation(param, Float64(factor), Float64(fraction), Float64(t_start), Float64(t_end))

"""
    Perturbation(schedule=ConstantDose(0.0); effects=[], gene_perturbations=[])

A dose schedule with the drug effects it drives, plus optional genetic
perturbations.
"""
struct Perturbation
    schedule::DoseSchedule
    effects::Vector{DrugEffect}
    gene_perturbations::Vector{GenePerturbation}
end
Perturbation(schedule::DoseSchedule=ConstantDose(0.0); effects=DrugEffect[], gene_perturbations=GenePerturbation[]) =
    Perturbation(schedule, collect(DrugEffect, effects), collect(GenePerturbation, gene_perturbations))

# ---- compiled forms ---------------------------------------------------------
struct DeathHazardC
    h_max::Float64; EC50::Float64; m::Float64; protect::Int; K::Float64; q::Float64; conc::Bool
end
struct RateModulationC
    idx::Int
    f::Function
end
struct GenePerturbationC
    idx::Int; factor::Float64; fraction::Float64; t_start::Float64; t_end::Float64
end
struct CompiledPerturbation
    schedule::DoseSchedule
    deaths::Vector{DeathHazardC}
    growth::Vector{GrowthInhibition}
    rates::Vector{RateModulationC}
    genes::Vector{GenePerturbationC}
end

function compile(pt::Perturbation, m::ReactionModel)
    deaths = DeathHazardC[]; growth = GrowthInhibition[]; rates = RateModulationC[]
    for e in pt.effects
        if e isa DeathHazard
            push!(deaths, DeathHazardC(e.h_max, e.EC50, e.m, e.protect === nothing ? 0 : speciesindex(m, e.protect), e.K, e.q, e.concentration))
        elseif e isa GrowthInhibition
            push!(growth, e)
        elseif e isa RateModulation
            push!(rates, RateModulationC(paramindex(m, e.param), e.f))
        end
    end
    genes = [GenePerturbationC(paramindex(m, g.param), g.factor, g.fraction, g.t_start, g.t_end) for g in pt.gene_perturbations]
    CompiledPerturbation(pt.schedule, deaths, growth, rates, genes)
end

@inline function hazard(e::DeathHazardC, d::Float64, x::AbstractVector{Int}, V::Float64)
    d <= 0.0 && return 0.0
    dm = d^e.m
    h = e.h_max * dm / (e.EC50^e.m + dm)
    if e.protect > 0
        c = e.conc ? x[e.protect] / V : float(x[e.protect])
        Kq = e.K^e.q
        h *= Kq / (Kq + c^e.q)
    end
    h
end

"""
    apply_effects!(cell, compiled_perturbation, p, d, t) -> (growth_multiplier, death_hazard)

Reset the cell's parameter vector to `p` and apply rate modulations and gene
perturbations; return the growth multiplier and the total drug death hazard.
"""
function apply_effects!(c::Cell, cp::Union{Nothing,CompiledPerturbation}, p::Vector{Float64}, d::Float64, t::Float64)
    copyto!(c.p, p)
    gmult = 1.0
    h = 0.0
    cp === nothing && return gmult, h
    for r in cp.rates
        c.p[r.idx] *= Float64(r.f(d))
    end
    if c.perturbed
        for g in cp.genes
            (g.t_start <= t < g.t_end) && (c.p[g.idx] *= g.factor)
        end
    end
    for gi in cp.growth
        gmult *= 1.0 / (1.0 + (d / gi.IC50)^gi.m)
    end
    for dh in cp.deaths
        h += hazard(dh, d, c.x, c.V)
    end
    gmult, h
end
