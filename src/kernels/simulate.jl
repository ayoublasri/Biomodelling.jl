"""
    Trajectory

Result of [`simulate`](@ref): `t` (save times), `X` (time × species counts) and
`species` names. Index a species column with `traj[:mRNA]`.
"""
struct Trajectory
    t::Vector{Float64}
    X::Matrix{Int}
    species::Vector{Symbol}
end
Base.getindex(tr::Trajectory, s::Symbol) = tr.X[:, findfirst(==(s), tr.species)]
Base.show(io::IO, tr::Trajectory) = print(io, "Trajectory(", length(tr.t), " time points, species ", tr.species, ")")

"""
    simulate(model, x0, tspan; kernel=DirectSSA(), p=model.p0, V=1.0, copies=1,
             saveat=nothing, rng=Random.default_rng()) -> Trajectory

Simulate a single cell (fixed volume) and record the state at `saveat` times
(default: 200 equally spaced points).
"""
function simulate(m::ReactionModel, x0::AbstractVector{<:Integer}, tspan::Tuple{<:Real,<:Real};
                  kernel::AbstractKernel=DirectSSA(), p::AbstractVector{<:Real}=m.p0, V::Real=1.0,
                  copies::Integer=1, saveat=nothing, rng::AbstractRNG=Random.default_rng())
    t0, t1 = float(tspan[1]), float(tspan[2])
    ts = saveat === nothing ? collect(range(t0, t1; length=201)) : collect(Float64, saveat)
    x = Vector{Int}(x0)
    length(x) == nspecies(m) || throw(ArgumentError("x0 must have $(nspecies(m)) entries"))
    pp = Vector{Float64}(p)
    ws = Workspace(m)
    X = zeros(Int, length(ts), nspecies(m))
    t = t0
    for (i, ti) in enumerate(ts)
        ti > t && (t = simulate!(x, m, pp, Float64(V), Int(copies), t, ti, kernel, ws, rng))
        X[i, :] .= x
    end
    Trajectory(ts, X, speciesnames(m))
end

"""
    ensemble_final(model, x0, T, n; kernel, p, V, copies, rng, threads=true) -> Matrix{Int}

Final states (n × species) of `n` independent single-cell simulations of
duration `T`. Deterministic for a given `rng` regardless of threading.
"""
function ensemble_final(m::ReactionModel, x0::AbstractVector{<:Integer}, T::Real, n::Integer;
                        kernel::AbstractKernel=DirectSSA(), p::AbstractVector{<:Real}=m.p0, V::Real=1.0,
                        copies::Integer=1, rng::AbstractRNG=Random.default_rng(), threads::Bool=true)
    pp = Vector{Float64}(p)
    seeds = rand(rng, UInt64, n)
    out = zeros(Int, n, nspecies(m))
    nth = threads ? Threads.nthreads() : 1
    if nth > 1 && n > 1
        @sync for chunk in _chunks(n, nth)
            Threads.@spawn begin
                # `local` keeps these task-private: without it they would be captured from
                # the enclosing function scope and shared between tasks.
                local tws = Workspace(m)
                local tx = Vector{Int}(undef, nspecies(m))
                for i in chunk
                    copyto!(tx, x0)
                    simulate!(tx, m, pp, Float64(V), Int(copies), 0.0, float(T), kernel, tws, Xoshiro(seeds[i]))
                    out[i, :] .= tx
                end
            end
        end
    else
        ws = Workspace(m)
        x = Vector{Int}(undef, nspecies(m))
        for i in 1:n
            copyto!(x, x0)
            simulate!(x, m, pp, Float64(V), Int(copies), 0.0, float(T), kernel, ws, Xoshiro(seeds[i]))
            out[i, :] .= x
        end
    end
    out
end

"""
    _chunks(n, k) -> Vector{UnitRange{Int}}

Split `1:n` into at most `k` contiguous, non-empty ranges.
"""
function _chunks(n::Int, k::Int)
    k = max(1, min(k, n))
    [(round(Int, (j - 1) * n / k) + 1):round(Int, j * n / k) for j in 1:k]
end
