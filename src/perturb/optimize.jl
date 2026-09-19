# Treatment outcomes and schedule optimisation.

"""
    cumulative_dose(result::PopulationResult) -> Float64

Dose actually applied in a simulation, integrated over the recorded time points
(trapezoidal rule), including feedback schedules.
"""
function cumulative_dose(r::PopulationResult)
    acc = 0.0
    for i in 2:length(r.t)
        acc += 0.5 * (r.dose[i] + r.dose[i-1]) * (r.t[i] - r.t[i-1])
    end
    acc
end

_index_at(r::PopulationResult, t::Real) = clamp(searchsortedlast(r.t, t), 1, length(r.t))

"""
    net_growth_rate(result; from=result.t[1]) -> Float64

Net exponential growth rate of the population between `from` and the end of
the simulation, `log(N_end / N_from) / (t_end - from)`; an extinct population
returns `-Inf`.
"""
function net_growth_rate(r::PopulationResult; from::Real=r.t[1])
    i = _index_at(r, from)
    N1 = r.popsize[end]
    N1 <= 0 && return -Inf
    (log(N1) - log(r.popsize[i])) / (r.t[end] - r.t[i])
end

"""
    log_kill(result; from=result.t[1]) -> Float64

Maximal depth of response in log10 units: `log10(N_from / N_min)` with `N_min`
the smallest population recorded after `from` (`Inf` for extinction).
"""
function log_kill(r::PopulationResult; from::Real=r.t[1])
    i = _index_at(r, from)
    Nmin = minimum(view(r.popsize, i:length(r.popsize)))
    Nmin <= 0 ? Inf : log10(r.popsize[i] / Nmin)
end

"""
    time_to_progression(result; from=result.t[1], threshold=1.2, floor=0.0) -> (time, progressed)

Time at which the population first exceeds `threshold` times its running
minimum after `from` (a RECIST-like progression from the nadir) while being
at least `floor` times its size at `from`. Returns the end of the simulation
and `false` if progression is not reached (a censored outcome).
"""
function time_to_progression(r::PopulationResult; from::Real=r.t[1], threshold::Real=1.2, floor::Real=0.0)
    i0 = _index_at(r, from)
    N0 = r.popsize[i0]
    nadir = N0
    for i in i0:length(r.t)
        N = r.popsize[i]
        nadir = min(nadir, N)
        if i > i0 && N >= threshold * nadir && N >= floor * N0 && N > nadir
            return r.t[i], true
        end
    end
    r.t[end], false
end

"""
    extinction_probability(results) -> Float64

Fraction of simulations whose population is extinct at the end.
"""
extinction_probability(results) = count(r -> r.popsize[end] <= 0, results) / length(results)

"""
    ScheduleOptimum

Result of [`optimize_schedule`](@ref): `params` and `value` of the best schedule,
the `names` of the parameters, and `table`, a vector of `(params, value)` for
every evaluated candidate (grid and refinement).
"""
struct ScheduleOptimum
    params::Vector{Float64}
    value::Float64
    names::Vector{Symbol}
    table::Vector{Tuple{Vector{Float64},Float64}}
end
Base.show(io::IO, o::ScheduleOptimum) = print(io, "ScheduleOptimum(", join(("$n = $(round(v; sigdigits=4))" for (n, v) in zip(o.names, o.params)), ", "), "; value = ", round(o.value; sigdigits=4), ", ", length(o.table), " evaluations)")

_lhs(rng, n, d) = [(randperm(rng, n) .- rand(rng, n)) ./ n for _ in 1:d]   # one stratified column per dimension

"""
    optimize_schedule(f, lower, upper; names, n_grid=5, max_grid=256, seeds=1:1, refine=true,
                      maxiter=60, log_scale=true, rng=Random.default_rng(), verbose=false)

Minimise a treatment objective over schedule parameters in the box
`[lower, upper]`. `f(θ, seed)` must return the objective for the parameter
vector `θ` (for instance the net growth rate, the negative time to
progression, or a penalised combination with a cumulative-dose constraint);
it is averaged over `seeds`, which are passed to every candidate so that
schedules are compared with common random numbers. The search evaluates a
log-spaced (or linear) grid of `n_grid` points per dimension (a Latin
hypercube of `max_grid` points when the full grid is larger) and then refines
the best point with Nelder-Mead in a bounded transform of the box.
"""
function optimize_schedule(f, lower::AbstractVector, upper::AbstractVector; names=[Symbol("p", i) for i in eachindex(lower)],
                           n_grid::Int=5, max_grid::Int=256, seeds=1:1, refine::Bool=true, maxiter::Int=60,
                           log_scale=true, rng::AbstractRNG=Random.default_rng(), verbose::Bool=false)
    d = length(lower)
    lo = Float64.(lower); hi = Float64.(upper)
    logs = log_scale isa Bool ? fill(log_scale, d) : collect(Bool, log_scale)
    tolog(θ) = [logs[i] ? log(θ[i]) : θ[i] for i in 1:d]
    fromlog(z) = [logs[i] ? exp(z[i]) : z[i] for i in 1:d]
    zlo = tolog(lo); zhi = tolog(hi)
    table = Tuple{Vector{Float64},Float64}[]
    objective(θ) = begin
        v = mean(f(θ, s) for s in seeds)
        push!(table, (copy(θ), v))
        verbose && println("  ", join(("$n = $(round(x; sigdigits=4))" for (n, x) in zip(names, θ)), ", "), " → ", round(v; sigdigits=4))
        v
    end
    # grid or Latin hypercube
    if n_grid^d <= max_grid
        axes = [range(zlo[i], zhi[i]; length=n_grid) for i in 1:d]
        candidates = [collect(Float64, c) for c in Iterators.product(axes...)]
    else
        u = _lhs(rng, max_grid, d)
        candidates = [[zlo[i] + (zhi[i] - zlo[i]) * u[i][k] for i in 1:d] for k in 1:max_grid]
    end
    best = nothing; bestv = Inf
    for z in candidates
        v = objective(fromlog(z))
        if v < bestv
            bestv = v; best = z
        end
    end
    if refine && d >= 1
        # bounded transform z = zlo + (zhi - zlo) * logistic(u)
        logistic(u) = 1 / (1 + exp(-u))
        logit(p) = log(p / (1 - p))
        u0 = [logit(clamp((best[i] - zlo[i]) / (zhi[i] - zlo[i]), 0.02, 0.98)) for i in 1:d]
        g(u) = objective(fromlog([zlo[i] + (zhi[i] - zlo[i]) * logistic(u[i]) for i in 1:d]))
        res = Optim.optimize(g, u0, Optim.NelderMead(), Optim.Options(iterations=maxiter, g_tol=1e-6))
        u = Optim.minimizer(res)
        z = [zlo[i] + (zhi[i] - zlo[i]) * logistic(u[i]) for i in 1:d]
        v = Optim.minimum(res)
        if v < bestv
            bestv = v; best = z
        end
    end
    ScheduleOptimum(fromlog(best), bestv, collect(Symbol, names), table)
end
