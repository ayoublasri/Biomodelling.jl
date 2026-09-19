"""
    timelapse(result, species; interval=1.0, σ=0.0, n_lineages=10, normalize=:concentration, rng)

Time-lapse reporter observation: follow the ancestral lines of `n_lineages`
cells of the final record,
sample the species (concentration or count) every `interval` time units and add
Gaussian measurement noise of standard deviation `σ` (multiplicative, relative
to the signal, if `σ` is negative: `|σ|` gives the coefficient of variation).
"""
function timelapse(r::PopulationResult, species; interval::Real=1.0, σ::Real=0.0, n_lineages::Int=10,
                   normalize::Symbol=:concentration, rng::AbstractRNG=Random.default_rng())
    s = _sidx(r.model, species)
    finals = r.ids[end]
    out = NamedTuple[]
    every = max(1, round(Int, interval / (r.t[2] - r.t[1])))
    for f in finals[1:min(n_lineages, length(finals))]
        tr = follow_lineage(r, f)
        isempty(tr.t) && continue
        idx = 1:every:length(tr.t)
        y = normalize == :concentration ? tr.X[idx, s] ./ tr.V[idx] : float.(tr.X[idx, s])
        if σ > 0
            y = y .+ σ .* randn(rng, length(y))
        elseif σ < 0
            y = y .* (1 .+ abs(σ) .* randn(rng, length(y)))
        end
        push!(out, (t = tr.t[idx], y = y, V = tr.V[idx], ids = tr.ids[idx], founder = f))
    end
    out
end
