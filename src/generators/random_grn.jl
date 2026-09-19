# ---------------------------------------------------------------------------
# Random gene regulatory networks
# ---------------------------------------------------------------------------

"""
    random_grn(n_genes; n_activations=n_genes, n_inhibitions=n_genes÷2, topology=:er,
               telegraph=false, k_tx=10.0, k_dm=1.0, K=5.0, n=2.0, basal=0.05,
               k_on=0.5, k_off=0.5, self_regulation=false, rng) -> (model, adjacency)

Build a random regulatory network of `n_genes` genes. `adjacency[i, j] = +1`
means gene `i` activates gene `j`, `-1` that it inhibits it. Regulation acts on
the transcription rate (constitutive genes, `telegraph = false`) or on the
promoter activation rate (`telegraph = true`), through concentration-based Hill
functions combined with AND logic. Scalar rate arguments may be given as a
`(lo, hi)` tuple to draw per-gene values log-uniformly.
"""
function random_grn(n_genes::Int; n_activations::Int=n_genes, n_inhibitions::Int=n_genes ÷ 2, topology::Symbol=:er,
                    telegraph::Bool=false, k_tx=10.0, k_dm=1.0, K=5.0, n=2.0, basal=0.05, k_on=0.5, k_off=0.5,
                    self_regulation::Bool=false, rng::AbstractRNG=Random.default_rng())
    draw(v) = v isa Tuple ? exp(log(v[1]) + rand(rng) * (log(v[2]) - log(v[1]))) : Float64(v)
    adj = zeros(Int, n_genes, n_genes)
    nedges = n_activations + n_inhibitions
    maxedges = n_genes * (n_genes - (self_regulation ? 0 : 1))
    nedges <= maxedges || throw(ArgumentError("too many edges for $n_genes genes"))
    outdeg = zeros(Int, n_genes)
    signs = shuffle(rng, vcat(fill(1, n_activations), fill(-1, n_inhibitions)))
    for sgn in signs
        while true
            src = if topology == :scale_free
                w = outdeg .+ 1
                sample(rng, 1:n_genes, Weights(w))
            else
                rand(rng, 1:n_genes)
            end
            dst = rand(rng, 1:n_genes)
            (src == dst && !self_regulation) && continue
            adj[src, dst] != 0 && continue
            adj[src, dst] = sgn
            outdeg[src] += 1
            break
        end
    end
    rx = Reaction[]
    params = Dict{Symbol,Float64}()
    groups = Vector{Vector{Symbol}}()
    mr(j) = Symbol("mRNA_", j)
    for j in 1:n_genes
        acts = [mr(i) for i in 1:n_genes if adj[i, j] == 1]
        inhs = [mr(i) for i in 1:n_genes if adj[i, j] == -1]
        ktx = Symbol("k_tx_", j); kdm = Symbol("k_dm_", j)
        params[ktx] = draw(k_tx); params[kdm] = draw(k_dm)
        regulated = !(isempty(acts) && isempty(inhs))
        if telegraph
            Goff = Symbol("G", j, "_off"); Gon = Symbol("G", j, "_on")
            kon = Symbol("k_on_", j); koff = Symbol("k_off_", j)
            params[kon] = draw(k_on); params[koff] = draw(k_off)
            kin = regulated ? Hill(kon; activators = acts, inhibitors = inhs, K = draw(K), n = draw(n), basal = draw(basal)) : MassAction(kon)
            push!(rx, Reaction("activation_$j", [Goff], [Gon], kin))
            push!(rx, Reaction("inactivation_$j", [Gon], [Goff], MassAction(koff)))
            push!(rx, Reaction("transcription_$j", [Gon], [Gon, mr(j)], MassAction(ktx); volume = :proportional))
            push!(groups, [Goff, Gon])
        else
            kin = regulated ? Hill(ktx; activators = acts, inhibitors = inhs, K = draw(K), n = draw(n), basal = draw(basal)) : MassAction(ktx)
            push!(rx, Reaction("transcription_$j", [], [mr(j)], kin; volume = :proportional, copy_number = true))
        end
        push!(rx, Reaction("degradation_$j", [mr(j)], [], MassAction(kdm)))
    end
    species = telegraph ? reduce(vcat, [[Symbol("G", j, "_off"), Symbol("G", j, "_on"), mr(j)] for j in 1:n_genes]) : [mr(j) for j in 1:n_genes]
    model = ReactionModel(rx; params = params, species = species, promoters = groups)
    model, adj
end
