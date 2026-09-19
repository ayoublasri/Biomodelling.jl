abstract type GrowthModel end

"""
    ExponentialGrowth(rate; cv=0.0)

Exponential volume growth `dV/dt = λ V`. With `cv > 0` each newborn cell draws
its own rate from a log-normal distribution with mean `rate` and coefficient of
variation ≈ `cv` (a simple source of extrinsic noise).
"""
struct ExponentialGrowth <: GrowthModel
    rate::Float64
    cv::Float64
end
ExponentialGrowth(rate::Real; cv::Real=0.0) = ExponentialGrowth(Float64(rate), Float64(cv))

sample_growth_rate(g::ExponentialGrowth, rng::AbstractRNG) =
    g.cv > 0 ? g.rate * exp(g.cv * randn(rng) - g.cv^2 / 2) : g.rate
@inline grow(::ExponentialGrowth, V::Float64, λ::Float64, Δt::Float64) = V * exp(λ * Δt)
doubling_time(g::ExponentialGrowth) = log(2) / g.rate
