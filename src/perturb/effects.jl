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
    GrowthInhibition(; IC50, m=2.0, protect=nothing, K=1.0, q=2.0, concentration=true)

Multiplies the growth rate by `1 / (1 + (d / IC50)^m)`. With a protective
species (`protect`), the inhibition is reduced by the factor
`K^q / (K^q + c^q)` of that species' concentration (or count), so that cells
in the resistant state keep proliferating under drug.
"""
struct GrowthInhibition <: DrugEffect
    IC50::Float64
    m::Float64
    protect::Union{Nothing,Symbol}
    K::Float64
    q::Float64
    concentration::Bool
end
GrowthInhibition(; IC50, m=2.0, protect=nothing, K=1.0, q=2.0, concentration=true) =
    GrowthInhibition(Float64(IC50), Float64(m), protect, Float64(K), Float64(q), concentration)
GrowthInhibition(IC50::Real, m::Real) = GrowthInhibition(Float64(IC50), Float64(m), nothing, 1.0, 2.0, true)

"""
    GrowthCost(species; K, q=2.0, max_cost=0.5, concentration=true)

Fitness cost of a cell state, independent of the dose: the growth rate is
multiplied by `1 - max_cost · c^q / (K^q + c^q)`, where `c` is the concentration
(or count) of `species`, so that cells expressing a resistance protein grow up
to `max_cost` slower.
"""
struct GrowthCost <: DrugEffect
    species::Symbol
    K::Float64
    q::Float64
    max_cost::Float64
    concentration::Bool
end
GrowthCost(species::Symbol; K, q=2.0, max_cost=0.5, concentration=true) =
    GrowthCost(species, Float64(K), Float64(q), Float64(max_cost), concentration)

"""
    CycleSensitivity(; baseline=0.0, center=0.5, width=0.15)

Cell-cycle dependence of the drug death hazard. The hazard of every
[`DeathHazard`](@ref) is multiplied by

    baseline + (1 - baseline) exp(-((φ - center)/width)^2 / 2),

where `φ` is the cell's progress through its division cycle (0 at birth, 1 at
division). The peak hazard is `h_max` at `φ = center` and falls to
`baseline · h_max` away from it, so `baseline = 1` recovers a cycle-independent
hazard and `baseline = 0` confines killing to a window around `center`. Setting
`center` to the replication set point represents an agent whose lesions are
converted into death during replication, such as a platinum drug.
"""
struct CycleSensitivity <: DrugEffect
    baseline::Float64
    center::Float64
    width::Float64
end
CycleSensitivity(; baseline=0.0, center=0.5, width=0.15) =
    CycleSensitivity(Float64(baseline), Float64(center), Float64(width))

"""
    SuicideConsumption(species; k, K_m=1.0)

Stoichiometric consumption of a protective protein by the drug, as for a
suicide enzyme such as MGMT, which is inactivated by the very lesion it
repairs. Lesions form at a rate proportional to the dose, so molecules of
`species` are removed at rate

    k · d · V · c / (c + K_m),

with dose `d`, volume `V` and protein concentration `c`. While the protein is
abundant the rate is set by lesion formation alone, so depletion follows the
*cumulative* exposure rather than the peak concentration; once the pool runs
low the rate falls with what remains. Resynthesis is whatever the reaction
model provides. This differs from scaling a first-order degradation rate by the
dose, which makes depletion track the peak concentration instead.
"""
struct SuicideConsumption <: DrugEffect
    species::Symbol
    k::Float64
    K_m::Float64
end
SuicideConsumption(species::Symbol; k, K_m=1.0) =
    SuicideConsumption(species, Float64(k), Float64(K_m))

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
struct GrowthCostC
    idx::Int; K::Float64; q::Float64; max_cost::Float64; conc::Bool
end
struct GenePerturbationC
    idx::Int; factor::Float64; fraction::Float64; t_start::Float64; t_end::Float64
end
struct GrowthInhibitionC
    IC50::Float64; m::Float64; protect::Int; K::Float64; q::Float64; conc::Bool
end
struct CycleSensitivityC
    baseline::Float64; center::Float64; width::Float64
end
struct SuicideConsumptionC
    idx::Int; k::Float64; K_m::Float64
end
struct CompiledPerturbation
    schedule::DoseSchedule
    deaths::Vector{DeathHazardC}
    growth::Vector{GrowthInhibitionC}
    rates::Vector{RateModulationC}
    genes::Vector{GenePerturbationC}
    costs::Vector{GrowthCostC}
    cycle::Vector{CycleSensitivityC}
    consume::Vector{SuicideConsumptionC}
end

function compile(pt::Perturbation, m::ReactionModel)
    deaths = DeathHazardC[]; growth = GrowthInhibitionC[]; rates = RateModulationC[]; costs = GrowthCostC[]
    cycle = CycleSensitivityC[]; consume = SuicideConsumptionC[]
    for e in pt.effects
        if e isa CycleSensitivity
            push!(cycle, CycleSensitivityC(e.baseline, e.center, e.width))
        elseif e isa SuicideConsumption
            push!(consume, SuicideConsumptionC(speciesindex(m, e.species), e.k, e.K_m))
        elseif e isa DeathHazard
            push!(deaths, DeathHazardC(e.h_max, e.EC50, e.m, e.protect === nothing ? 0 : speciesindex(m, e.protect), e.K, e.q, e.concentration))
        elseif e isa GrowthInhibition
            push!(growth, GrowthInhibitionC(e.IC50, e.m, e.protect === nothing ? 0 : speciesindex(m, e.protect), e.K, e.q, e.concentration))
        elseif e isa RateModulation
            push!(rates, RateModulationC(paramindex(m, e.param), e.f))
        elseif e isa GrowthCost
            push!(costs, GrowthCostC(speciesindex(m, e.species), e.K, e.q, e.max_cost, e.concentration))
        end
    end
    genes = [GenePerturbationC(paramindex(m, g.param), g.factor, g.fraction, g.t_start, g.t_end) for g in pt.gene_perturbations]
    CompiledPerturbation(pt.schedule, deaths, growth, rates, genes, costs, cycle, consume)
end

"""
    cycle_multiplier(cp, phi) -> Float64

Factor by which the drug death hazard is scaled for a cell at cycle progress
`phi`; `1.0` when no [`CycleSensitivity`](@ref) is in force.
"""
@inline function cycle_multiplier(cp::CompiledPerturbation, phi::Float64)
    isempty(cp.cycle) && return 1.0
    mult = 1.0
    for cs in cp.cycle
        z = (phi - cs.center) / cs.width
        mult *= cs.baseline + (1.0 - cs.baseline) * exp(-0.5 * z * z)
    end
    mult
end

"""
    apply_consumption!(cell, compiled_perturbation, d, dt)

Remove molecules of every species under [`SuicideConsumption`](@ref), one per
lesion repaired over the step, by a Poisson draw at the current dose.
"""
function apply_consumption!(c::Cell, cp::Union{Nothing,CompiledPerturbation}, d::Float64, dt::Float64)
    (cp === nothing || d <= 0.0 || dt <= 0.0) && return nothing
    for sc in cp.consume
        M = c.x[sc.idx]
        M <= 0 && continue
        conc = M / c.V
        rate = sc.k * d * c.V * conc / (conc + sc.K_m)
        rate > 0.0 || continue
        n = pois_rand(c.rng, rate * dt)
        c.x[sc.idx] = max(0, M - n)
    end
    nothing
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
    apply_effects!(cell, compiled_perturbation, p, d, t, phi=0.5) -> (growth_multiplier, death_hazard)

Reset the cell's parameter vector to `p` and apply rate modulations and gene
perturbations; return the growth multiplier and the total drug death hazard.
`phi` is the cell's progress through its division cycle, used by
[`CycleSensitivity`](@ref).
"""
function apply_effects!(c::Cell, cp::Union{Nothing,CompiledPerturbation}, p::Vector{Float64}, d::Float64, t::Float64, phi::Float64=0.5)
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
        inhib = 1.0 - 1.0 / (1.0 + (d / gi.IC50)^gi.m)
        if gi.protect > 0
            conc = gi.conc ? c.x[gi.protect] / c.V : float(c.x[gi.protect])
            Kq = gi.K^gi.q
            inhib *= Kq / (Kq + conc^gi.q)
        end
        gmult *= 1.0 - inhib
    end
    for gc in cp.costs
        conc = gc.conc ? c.x[gc.idx] / c.V : float(c.x[gc.idx])
        cq = conc^gc.q
        gmult *= 1.0 - gc.max_cost * cq / (gc.K^gc.q + cq)
    end
    for dh in cp.deaths
        h += hazard(dh, d, c.x, c.V)
    end
    h > 0.0 && (h *= cycle_multiplier(cp, phi))
    gmult, h
end
