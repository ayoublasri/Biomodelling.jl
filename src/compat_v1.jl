# ---------------------------------------------------------------------------
# Compatibility layer for the v0.3 (2019-2022) API. Deprecated: will be removed in v3.
# ---------------------------------------------------------------------------

"""
    Donne(model, initiale, T, tau, NoC, growth_data, epsilon=0.03; trans_index=Int[])

Deprecated v1 model container. `model` is a tuple of NamedTuple reactions with
fields `name, rate, reactants, products, coeff_rea, coeff_pro`; `initiale` is a
two-column matrix of species symbols and counts (a `:NULL` row is allowed).
Reactions whose name contains `"act"`, `"inhib"`, `"comb_a"` or `"comb_i"` use
the original Hill-type propensities with `rate = [k, n, K]`. Reactions listed in
`trans_index` scale with cell volume. Prefer [`ReactionModel`](@ref).
"""
struct Donne
    M::Int
    N::Int
    T::Float64
    tau::Float64
    NoC::Int
    growth_rate::Float64
    epsilon::Float64
    species::Vector{Symbol}
    X::Vector{Int}
    model::Any
    rmodel::ReactionModel
    x0::Vector{Int}
    trans_index::Vector{Int}
end

function _v1_convert(model, trans_index::Vector{Int})
    rxns = Reaction[]
    params = Dict{Symbol,Float64}()
    for (i, r) in enumerate(model)
        name = String(r.name)
        lname = lowercase(name)
        reactants = collect(Symbol, r.reactants); products = collect(Symbol, r.products)
        cr = collect(Int, r.coeff_rea); cp = collect(Int, r.coeff_pro)
        volume = i in trans_index ? :proportional : :none
        regulated = occursin("act", lname) || occursin("inhib", lname) || occursin("comb_a", lname) || occursin("comb_i", lname)
        if regulated
            net = Dict{Symbol,Int}()
            for (s, c) in zip(reactants, cr)
                s in _NULLS && continue
                net[s] = get(net, s, 0) - c
            end
            for (s, c) in zip(products, cp)
                s in _NULLS && continue
                net[s] = get(net, s, 0) + c
            end
            regs = [s for s in reactants if !(s in _NULLS)]
            comb = occursin("comb", lname)
            comb || (regs = regs[1:1])
            mode = (occursin("inhib", lname) || occursin("comb_i", lname)) ? :inhibit : :activate
            rv = collect(Float64, r.rate isa Number ? [r.rate] : r.rate)
            k = rv[1]; n = length(rv) >= 2 ? rv[2] : 1.0; K = length(rv) >= 3 ? rv[3] : 1.0
            vf = (length(rv) >= 4 && volume == :none) ? rv[4] : 1.0
            tag = _sanitize(name) * "_$i"
            pk = Symbol("k_", tag); pn = Symbol("n_", tag); pK = Symbol("K_", tag)
            params[pk] = k; params[pn] = n; params[pK] = K
            regnames = copy(regs)
            kin = Hill(pk; activators = mode == :activate ? regnames : Symbol[],
                       inhibitors = mode == :inhibit ? regnames : Symbol[], K = pK, n = pn, basal = 0.0)
            prods = [s => c for (s, c) in net if c != 0]
            # a regulated reaction has no mass-action reactants; net changes are expressed as products
            push!(rxns, Reaction(name, Symbol[], prods, kin; volume = volume))
            vf == 1.0 || (params[pk] = k * vf)
        else
            pk = Symbol("k_", _sanitize(name), "_$i")
            params[pk] = Float64(r.rate isa Number ? r.rate : r.rate[1])
            rea = [s => c for (s, c) in zip(reactants, cr) if !(s in _NULLS)]
            pro = [s => c for (s, c) in zip(products, cp) if !(s in _NULLS)]
            push!(rxns, Reaction(name, rea, pro, MassAction(pk); volume = volume))
        end
    end
    order = Symbol[]
    for r in model, s in vcat(collect(Symbol, r.reactants), collect(Symbol, r.products))
        (s in _NULLS || s in order) || push!(order, s)
    end
    # promoter groups: species whose names contain "on"/"off" (v1 convention), paired by order
    ons = [s for s in order if occursin("on", lowercase(String(s))) && !occursin("off", lowercase(String(s)))]
    offs = [s for s in order if occursin("off", lowercase(String(s)))]
    groups = Vector{Vector{Symbol}}()
    for (a, b) in zip(ons, offs)
        push!(groups, [b, a])
    end
    ReactionModel(rxns; params = params, species = order, promoters = groups)
end

function Donne(model, initiale, T, tau, NoC, growth_data, epsilon=0.03; trans_index=Int[])
    Base.depwarn("`Donne` is the deprecated v1 API; build a `ReactionModel` instead", :Donne)
    rmodel = _v1_convert(model, collect(Int, trans_index))
    species = vcat([:NULL], speciesnames(rmodel))
    x0 = zeros(Int, nspecies(rmodel))
    for i in 1:size(initiale, 1)
        s = initiale[i, 1]
        s in _NULLS && continue
        x0[speciesindex(rmodel, Symbol(s))] = Int(initiale[i, 2])
    end
    growth_rate = growth_data isa Real ? Float64(growth_data) : growth_estimate(growth_data)
    Donne(nreactions(rmodel), length(species), Float64(T), Float64(tau), Int(NoC), growth_rate, Float64(epsilon),
          species, vcat([0], x0), model, rmodel, x0, collect(Int, trans_index))
end
Donne(model, matrix::AbstractMatrix{<:Number}, initiale::AbstractMatrix, T, tau, NoC, growth_data, epsilon=0.03) =
    error("the v1 `comb_rea` matrix form is not supported by the compatibility layer; use `Hill` kinetics")

"""
    growth_estimate(growth_data) -> Float64

Exponential growth rate fitted by log-linear least squares to columns 1 (time)
and 5 (count) of a matrix, as in v1.
"""
function growth_estimate(growth_data::AbstractMatrix)
    t = Float64.(growth_data[:, 1]); y = Float64.(growth_data[:, 5])
    keep = y .> 0
    t = t[keep]; ly = log.(y[keep])
    tm = mean(t); lm = mean(ly)
    sum((t .- tm) .* (ly .- lm)) / sum((t .- tm) .^ 2)
end

function _v1_single(data::Donne, kernel::AbstractKernel; rng::AbstractRNG=Random.default_rng())
    ts = collect(0.0:data.tau:data.T)
    tr = simulate(data.rmodel, data.x0, (0.0, data.T); kernel, saveat = ts, rng)
    X = hcat(zeros(Int, length(ts)), tr.X)
    ts, X
end

"""
    ssa(data::Donne) -> (t, X)    (deprecated)
"""
ssa(data::Donne; rng::AbstractRNG=Random.default_rng()) = _v1_single(data, DirectSSA(); rng)
"""
    tauleap(data::Donne) -> (t, X)    (deprecated)
"""
tauleap(data::Donne; rng::AbstractRNG=Random.default_rng()) = _v1_single(data, TauLeap(data.tau); rng)
"""
    tauleapswitch(data::Donne, ssa_steps=10) -> (t, X)    (deprecated; hybrid tau-leap/SSA)
"""
tauleapswitch(data::Donne, ssa_steps::Integer=10; rng::AbstractRNG=Random.default_rng()) = _v1_single(data, HybridSSATau(data.tau); rng)
"""
    adaptive_tauleap(data::Donne) -> (t, X)    (deprecated)
"""
adaptive_tauleap(data::Donne; rng::AbstractRNG=Random.default_rng()) = _v1_single(data, AdaptiveTauLeap(ε = data.epsilon); rng)

_v1_kernel(alg, data::Donne) = alg === ssa ? DirectSSA() :
                               alg === tauleap ? TauLeap(data.tau) :
                               alg === adaptive_tauleap ? AdaptiveTauLeap(ε = data.epsilon) :
                               HybridSSATau(data.tau)

"""
    exponential_growth(data::Donne, div_noise, alg, Ni; rng) -> (t, V, X)
    exponential_growth(data::Donne, trans_index, div_noise, alg, Ni; rng)

Deprecated v1 population simulation: constant population of `data.NoC` cells,
sizer division at volume 2, binomial partitioning with division noise
`div_noise`, initial volumes `1 + Ni * rand()`. Returns time, a (time × cells)
volume matrix and a (time × cells × species) count array whose first species
column is `:NULL`.
"""
function exponential_growth(data::Donne, div_noise::Real, alg, Ni::Real; rng::AbstractRNG=Random.default_rng())
    settings = PopulationSettings(dt = data.tau, kernel = _v1_kernel(alg, data),
                                  growth = ExponentialGrowth(data.growth_rate), size_control = Sizer(2.0; cv = 0.0005),
                                  partitioning = BinomialPartition(σ = Float64(div_noise)), control = ConstantN(),
                                  randomize_promoters = true, track_lineage = false)
    V0 = 1.0 .+ Ni .* rand(rng, data.NoC)
    res = simulate_population(data.rmodel, data.x0, data.NoC, (0.0, data.T); settings, rng, V0)
    nt = length(res.t)
    V = zeros(nt, data.NoC)
    X = zeros(Int, nt, data.NoC, data.N)
    for i in 1:nt
        n = min(size(res.counts[i], 1), data.NoC)
        V[i, 1:n] = res.volume[i][1:n]
        X[i, 1:n, 2:end] = res.counts[i][1:n, :]
    end
    res.t, V, X
end
function exponential_growth(data::Donne, trans_index::AbstractVector{<:Integer}, div_noise::Real, alg, Ni::Real; kwargs...)
    d = Donne(data.model, hcat(data.species, data.X), data.T, data.tau, data.NoC, data.growth_rate, data.epsilon;
              trans_index = collect(Int, trans_index))
    exponential_growth(d, div_noise, alg, Ni; kwargs...)
end
