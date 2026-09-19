"""
    DoseSchedule

Time course of an external perturbation (drug dose). Evaluate with
[`dose`](@ref)`(schedule, t)`. Concrete types: [`ConstantDose`](@ref),
[`PulsedDose`](@ref), [`PiecewiseDose`](@ref), [`BolusPK`](@ref), or any
function `t -> d` wrapped in `FunctionDose`.
"""
abstract type DoseSchedule end

struct ConstantDose <: DoseSchedule
    d::Float64
end

"""
    PulsedDose(d; on, off, start=0.0)

Dose `d` for `on` time units, then `0` for `off` time units, repeated from
`start` (drug holidays).
"""
struct PulsedDose <: DoseSchedule
    d::Float64
    on::Float64
    off::Float64
    start::Float64
end
PulsedDose(d::Real; on::Real, off::Real, start::Real=0.0) = PulsedDose(Float64(d), Float64(on), Float64(off), Float64(start))

"""
    PiecewiseDose(times, doses)

`doses[i]` applies on `[times[i], times[i+1])`; the dose is 0 before `times[1]`.
"""
struct PiecewiseDose <: DoseSchedule
    times::Vector{Float64}
    doses::Vector{Float64}
    function PiecewiseDose(times, doses)
        length(times) == length(doses) || throw(ArgumentError("times and doses must have equal length"))
        issorted(times) || throw(ArgumentError("times must be sorted"))
        new(Vector{Float64}(times), Vector{Float64}(doses))
    end
end

"""
    BolusPK(times, amounts, k_e)

One-compartment pharmacokinetics: each bolus of `amounts[i]` given at `times[i]`
decays exponentially with elimination rate `k_e`.
"""
struct BolusPK <: DoseSchedule
    times::Vector{Float64}
    amounts::Vector{Float64}
    k_e::Float64
end

struct FunctionDose <: DoseSchedule
    f::Function
end

dose(s::ConstantDose, t::Real) = s.d
function dose(s::PulsedDose, t::Real)
    t < s.start && return 0.0
    τ = mod(t - s.start, s.on + s.off)
    τ < s.on ? s.d : 0.0
end
function dose(s::PiecewiseDose, t::Real)
    i = searchsortedlast(s.times, t)
    i == 0 ? 0.0 : s.doses[i]
end
function dose(s::BolusPK, t::Real)
    c = 0.0
    for (ti, A) in zip(s.times, s.amounts)
        t >= ti && (c += A * exp(-s.k_e * (t - ti)))
    end
    c
end
dose(s::FunctionDose, t::Real) = Float64(s.f(t))
dose(::Nothing, t::Real) = 0.0
