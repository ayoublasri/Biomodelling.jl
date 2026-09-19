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
    PulsedDose(d; on, off, start=0.0, cycles=0)

Dose `d` for `on` time units, then `0` for `off` time units, repeated from
`start` (drug holidays). `cycles > 0` limits the number of on/off cycles, after
which the dose is 0 (a clinical regimen of `cycles` courses).
"""
struct PulsedDose <: DoseSchedule
    d::Float64
    on::Float64
    off::Float64
    start::Float64
    cycles::Int
end
PulsedDose(d::Real; on::Real, off::Real, start::Real=0.0, cycles::Integer=0) = PulsedDose(Float64(d), Float64(on), Float64(off), Float64(start), Int(cycles))
PulsedDose(d::Real, on::Real, off::Real, start::Real) = PulsedDose(Float64(d), Float64(on), Float64(off), Float64(start), 0)

"""
    AdaptiveDose(d; on_above=1.0, off_below=0.5, start=0.0, initial=true)

Feedback ("adaptive therapy") schedule in the spirit of Zhang et al. (2017):
the dose `d` is applied while the population is above `off_below` times its
size at `start` and re-applied once it has regrown above `on_above` times that
size. The schedule is stateful and is evaluated with the population size by
[`dose`](@ref)`(schedule, t, N)`, as the population loop does; `initial` sets
whether treatment starts at `start`.
"""
mutable struct AdaptiveDose <: DoseSchedule
    d::Float64
    on_above::Float64
    off_below::Float64
    start::Float64
    initial::Bool
    N_ref::Float64
    active::Bool
    t_last::Float64
end
AdaptiveDose(d::Real; on_above::Real=1.0, off_below::Real=0.5, start::Real=0.0, initial::Bool=true) =
    AdaptiveDose(Float64(d), Float64(on_above), Float64(off_below), Float64(start), initial, NaN, false, -Inf)

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
    s.cycles > 0 && t >= s.start + s.cycles * (s.on + s.off) && return 0.0
    τ = mod(t - s.start, s.on + s.off)
    τ < s.on ? s.d : 0.0
end
"""
    dose(schedule, t) -> Float64
    dose(schedule, t, N) -> Float64

Dose at time `t`; the three-argument form also receives the current population
size `N`, which feedback schedules ([`AdaptiveDose`](@ref)) use and all other
schedules ignore.
"""
dose(s::DoseSchedule, t::Real, N::Real) = dose(s, t)
function dose(s::AdaptiveDose, t::Real, N::Real)
    if t < s.start
        s.N_ref = NaN; s.active = false; s.t_last = t
        return 0.0
    end
    if isnan(s.N_ref) || t < s.t_last          # first evaluation, or a new simulation reusing the object
        s.N_ref = Float64(N); s.active = s.initial
    end
    s.t_last = t
    if s.active
        N <= s.off_below * s.N_ref && (s.active = false)
    else
        N >= s.on_above * s.N_ref && (s.active = true)
    end
    s.active ? s.d : 0.0
end
dose(s::AdaptiveDose, t::Real) = (t >= s.start && s.active) ? s.d : 0.0
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

"""
    daily_boluses(amount, days, k_e; start=0.0, interval=24.0) -> BolusPK

One-compartment pharmacokinetics of repeated boluses: `amount` given on each of
`days` (0-based day indices, e.g. `0:4` for days 1-5 of a cycle, or
`vcat(0:4, 28:32)` for two 28-day cycles), `interval` time units apart, from
`start`, each eliminated with rate `k_e` (`log(2) / half_life`).
"""
daily_boluses(amount::Real, days, k_e::Real; start::Real=0.0, interval::Real=24.0) =
    BolusPK([start + interval * d for d in days], fill(Float64(amount), length(days)), Float64(k_e))

"""
    cycle_days(days_on, cycle_length, n_cycles) -> Vector{Int}

Day indices (0-based) of a regimen of `n_cycles` cycles of `cycle_length` days
with the drug given on the first `days_on` days of each cycle.
"""
cycle_days(days_on::Integer, cycle_length::Integer, n_cycles::Integer) =
    [c * cycle_length + d for c in 0:n_cycles-1 for d in 0:days_on-1]

"""
    cumulative_dose(schedule, t0, t1; dt=0.05) -> Float64

Integral of the dose over `[t0, t1]` (for stateless schedules); use
`cumulative_dose(result)` for the dose actually applied in a simulation.
"""
function cumulative_dose(s::DoseSchedule, t0::Real, t1::Real; dt::Real=0.05)
    acc = 0.0
    t = float(t0)
    while t < t1
        acc += dose(s, t) * min(dt, t1 - t)
        t += dt
    end
    acc
end
