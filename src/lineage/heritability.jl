# ---------------------------------------------------------------------------
# Heritability and memory statistics from the lineage table
# ---------------------------------------------------------------------------

_sidx(m::ReactionModel, s) = s isa Integer ? Int(s) : speciesindex(m, s)

function _end_value(lt::LineageTable, id::Int, s::Int, normalize::Symbol)
    x = lt.x_end[id]
    isempty(x) && return NaN
    normalize == :concentration ? x[s] / lt.V_end[id] : float(x[s])
end

function _pair_values(r::PopulationResult, s::Int, pairs, normalize::Symbol; divided_only::Bool=true)
    lt = r.lineage
    a = Float64[]; b = Float64[]
    for (i, j) in pairs
        if divided_only
            (lt.fate[i] == :divided && lt.fate[j] == :divided) || continue
        else
            (isempty(lt.x_end[i]) || isempty(lt.x_end[j])) && continue
        end
        push!(a, _end_value(lt, i, s, normalize)); push!(b, _end_value(lt, j, s, normalize))
    end
    a, b
end

"""
    heritability(result, species; relation=:mother_daughter, normalize=:concentration,
                 divided_only=true) -> (r, n)

Pearson correlation of expression at division between related cells:
`relation` is `:mother_daughter`, `:sisters` or `:cousins`. With `divided_only`
only cells whose life ended by division are used, so both members of a pair are
compared at the same cell-cycle stage.
"""
function heritability(r::PopulationResult, species; relation::Symbol=:mother_daughter,
                      normalize::Symbol=:concentration, divided_only::Bool=true)
    s = _sidx(r.model, species)
    pairs = relation == :mother_daughter ? mother_daughter_pairs(r.lineage) :
            relation == :sisters ? sister_pairs(r.lineage) :
            relation == :cousins ? cousin_pairs(r.lineage) :
            throw(ArgumentError("relation must be :mother_daughter, :sisters or :cousins"))
    a, b = _pair_values(r, s, pairs, normalize; divided_only)
    n = length(a)
    n < 3 && return (r = NaN, n = n)
    (std(a) == 0 || std(b) == 0) && return (r = NaN, n = n)
    (r = cor(a, b), n = n)
end

"""
    lineage_autocorrelation(result, species; max_generations=6, normalize=:concentration) -> Vector{Float64}

Correlation of expression at division between a cell and its ancestor `g`
generations earlier, for `g = 1:max_generations`.
"""
function lineage_autocorrelation(r::PopulationResult, species; max_generations::Int=6, normalize::Symbol=:concentration)
    s = _sidx(r.model, species)
    lt = r.lineage
    out = fill(NaN, max_generations)
    divided = findall(==(:divided), lt.fate)
    for g in 1:max_generations
        a = Float64[]; b = Float64[]
        for id in divided
            anc = id
            ok = true
            for _ in 1:g
                anc = lt.parent[anc]
                if anc == 0 || lt.fate[anc] != :divided
                    ok = false; break
                end
            end
            ok || continue
            push!(a, _end_value(lt, id, s, normalize)); push!(b, _end_value(lt, anc, s, normalize))
        end
        if length(a) >= 3 && std(a) > 0 && std(b) > 0
            out[g] = cor(a, b)
        end
    end
    out
end

"""
    memory_timescale(result, species; max_generations=6) -> (generations, time, correlations)

Fit `r_g = exp(-g / τ)` to the lineage autocorrelation and return the memory
timescale `τ` in generations and in time units (using the mean cell-cycle time
observed in the lineage table).
"""
function memory_timescale(r::PopulationResult, species; max_generations::Int=6, normalize::Symbol=:concentration)
    rg = lineage_autocorrelation(r, species; max_generations, normalize)
    gs = Float64[]; ys = Float64[]
    for (g, v) in enumerate(rg)
        (isfinite(v) && v > 0) && (push!(gs, g); push!(ys, log(v)))
    end
    lt = r.lineage
    cyc = [lt.end_time[i] - lt.birth_time[i] for i in eachindex(lt.id) if lt.fate[i] == :divided && lt.generation[i] > 0]
    Tc = isempty(cyc) ? NaN : mean(cyc)
    if length(gs) < 2
        τ = length(gs) == 1 ? -gs[1] / ys[1] : NaN
    else
        # least squares through the origin: log r = -g/τ
        τ = -sum(gs .* gs) / sum(gs .* ys)
    end
    (generations = τ, time = τ * Tc, correlations = rg, cycle_time = Tc)
end

"""
    follow_lineage(result, id) -> (t, X, V, ids)

Reconstruct the single-cell trace of the ancestral line of cell `id` (normally a
cell present in the last record): at every recorded time the ancestor alive at
that time is reported, so the trace runs from the start of the simulation to the
end of the cell's life.
"""
function follow_lineage(r::PopulationResult, id::Integer)
    lt = r.lineage
    chain = reverse(lineage_of(lt, Int(id)))          # founder first
    births = [lt.birth_time[c] for c in chain]
    ts = Float64[]; X = Vector{Vector{Int}}(); V = Float64[]; ids = Int[]
    for i in eachindex(r.t)
        k = searchsortedlast(births, r.t[i])
        k == 0 && continue
        cur = chain[k]
        row = findfirst(==(cur), r.ids[i])
        row === nothing && (k == length(chain) ? break : continue)
        push!(ts, r.t[i]); push!(X, r.counts[i][row, :]); push!(V, r.volume[i][row]); push!(ids, cur)
    end
    (t = ts, X = isempty(X) ? zeros(Int, 0, nspecies(r.model)) : permutedims(reduce(hcat, X)), V = V, ids = ids)
end

"""
    noise_decomposition(result, species; n_lineages=20, burnin=0.2) -> NamedTuple

Coefficient of variation squared of a species' concentration in the population
perspective (across cells at the final record) and in the single-lineage
perspective (along time in the ancestral lines of `n_lineages` cells of the
final record, after discarding the first `burnin` fraction of the records).
"""
function noise_decomposition(r::PopulationResult, species; n_lineages::Int=20, burnin::Float64=0.2, normalize::Symbol=:concentration)
    s = _sidx(r.model, species)
    i0 = max(1, round(Int, burnin * length(r.t)))
    pop = normalize == :concentration ? r.counts[end][:, s] ./ r.volume[end] : float.(r.counts[end][:, s])
    cv2_pop = var(pop) / mean(pop)^2
    finals = r.ids[end]
    cv2s = Float64[]
    for f in finals[1:min(n_lineages, length(finals))]
        tr = follow_lineage(r, f)
        length(tr.t) < 10 && continue
        y = normalize == :concentration ? tr.X[i0:end, s] ./ tr.V[i0:end] : float.(tr.X[i0:end, s])
        (length(y) < 5 || mean(y) <= 0) && continue          # CV² undefined for silent lineages
        push!(cv2s, var(y) / mean(y)^2)
    end
    (population = cv2_pop, lineage = isempty(cv2s) ? NaN : mean(cv2s), lineage_values = cv2s, n_lineages = length(cv2s))
end
