abstract type SizeControl end

"""
    Sizer(V_div; cv=0.0)

Divide when the volume reaches `V_div` (times a per-cell log-normal-like noise
factor of coefficient of variation `cv`, drawn at birth).
"""
struct Sizer <: SizeControl
    V_div::Float64
    cv::Float64
end
Sizer(V_div::Real; cv::Real=0.0) = Sizer(Float64(V_div), Float64(cv))

"""
    Adder(Δ; cv=0.0)

Divide after adding volume `Δ` since birth (adder size control).
"""
struct Adder <: SizeControl
    Δ::Float64
    cv::Float64
end
Adder(Δ::Real; cv::Real=0.0) = Adder(Float64(Δ), Float64(cv))

"""
    AgeTimer(T; cv=0.0)

Divide at age `T` (times noise), independent of size.
"""
struct AgeTimer <: SizeControl
    T::Float64
    cv::Float64
end
AgeTimer(T::Real; cv::Real=0.0) = AgeTimer(Float64(T), Float64(cv))

"""
    Replication(fraction)

Replicate genes (double promoter counts and the `copies` multiplier) when the
cell has completed `fraction` of its cycle (0 < fraction < 1).
"""
struct Replication
    fraction::Float64
    function Replication(f::Real)
        0 < f < 1 || throw(ArgumentError("replication fraction must be in (0, 1)"))
        new(Float64(f))
    end
end

_noisy(v::Float64, cv::Float64, rng::AbstractRNG) = cv > 0 ? v * max(0.05, 1.0 + cv * randn(rng)) : v

sample_target(sc::Sizer, Vb::Float64, rng::AbstractRNG) = _noisy(sc.V_div, sc.cv, rng)
sample_target(sc::Adder, Vb::Float64, rng::AbstractRNG) = Vb + _noisy(sc.Δ, sc.cv, rng)
sample_target(sc::AgeTimer, Vb::Float64, rng::AbstractRNG) = _noisy(sc.T, sc.cv, rng)

@inline should_divide(::Union{Sizer,Adder}, c::Cell, t::Float64) = c.V >= c.div_target
@inline should_divide(::AgeTimer, c::Cell, t::Float64) = (t - c.birth_time) >= c.div_target

@inline cycle_progress(::Union{Sizer,Adder}, c::Cell, t::Float64) =
    (c.V - c.V_birth) / max(c.div_target - c.V_birth, 1e-12)
@inline cycle_progress(::AgeTimer, c::Cell, t::Float64) = (t - c.birth_time) / c.div_target

typical_birth_volume(sc::Sizer) = sc.V_div / 2
typical_birth_volume(sc::Adder) = sc.Δ
typical_birth_volume(::AgeTimer) = 1.0

function replicate!(c::Cell, m::ReactionModel)
    for g in m.promoter_groups, s in g
        c.x[s] *= 2
    end
    c.copies = 2
    c.replicated = true
    c
end
