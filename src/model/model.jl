# ---------------------------------------------------------------------------
# Compiled reaction model
# ---------------------------------------------------------------------------

"""
    Species(name; kind=:molecule)

A chemical species. `kind` is `:molecule` (counts are partitioned binomially at
division) or `:promoter` (a gene state, inherited by daughters; see the
`promoters` keyword of [`ReactionModel`](@ref)).
"""
struct Species
    name::Symbol
    kind::Symbol
end
Species(name::Symbol; kind::Symbol=:molecule) = Species(name, kind)

struct MassActionC
    rate::Int
    reactants::Vector{Tuple{Int,Int}}
    vol_exp::Float64
    copy_number::Bool
end

struct HillC
    rate::Int
    reactants::Vector{Tuple{Int,Int}}
    vol_exp::Float64
    copy_number::Bool
    regs::Vector{Int}
    activate::Vector{Bool}
    K::Vector{Int}
    n::Vector{Int}
    basal::Int
    and_logic::Bool
    conc::Bool
end

struct CustomC
    f::Function
    params::Vector{Int}
    reactants::Vector{Tuple{Int,Int}}
    vol_exp::Float64
    copy_number::Bool
end

const CompiledKinetics = Union{MassActionC,HillC,CustomC}

"""
    ReactionModel(reactions; params, species=nothing, promoters=[])

Compile a vector of [`Reaction`](@ref)s into a simulation-ready model.

* `params`: `NamedTuple` or `Dict` giving values for every named parameter.
* `species`: optional ordering of species names (species not listed are
  appended in order of appearance).
* `promoters`: vector of groups of species names that together represent the
  states of one gene's promoter(s), e.g. `[[:G_off, :G_on]]`. Promoter counts are
  doubled at gene replication and inherited (not binomially partitioned) at
  division.

Fields of interest: `species`, `reactions`, `pnames`, `p0` (default parameter
vector), `stoich` (per-reaction sparse net change), `depgraph`.
"""
struct ReactionModel
    species::Vector{Species}
    sidx::Dict{Symbol,Int}
    reactions::Vector{Reaction}
    kinetics::Vector{CompiledKinetics}
    stoich::Vector{Vector{Tuple{Int,Int}}}
    reactant_stoich::Vector{Vector{Tuple{Int,Int}}}
    depgraph::Vector{Vector{Int}}
    pnames::Vector{Symbol}
    pidx::Dict{Symbol,Int}
    p0::Vector{Float64}
    hor::Vector{Int}
    mrr::Vector{Int}
    promoter_groups::Vector{Vector{Int}}
    is_promoter::Vector{Bool}
end

function _sanitize(name::AbstractString)
    s = replace(name, r"[^A-Za-z0-9]+" => "_")
    isempty(s) ? "rxn" : s
end

function ReactionModel(reactions::AbstractVector{Reaction}; params=NamedTuple(),
                       species=nothing, promoters=Vector{Vector{Symbol}}())
    rxns = collect(Reaction, reactions)
    isempty(rxns) && throw(ArgumentError("a model needs at least one reaction"))
    # --- species -----------------------------------------------------------
    names = Symbol[]
    addname!(s::Symbol) = (s in names || push!(names, s); nothing)
    species === nothing || foreach(addname!, species)
    for r in rxns
        foreach(x -> addname!(first(x)), r.reactants)
        foreach(x -> addname!(first(x)), r.products)
        r.kinetics isa Hill && foreach(addname!, r.kinetics.regulators)
    end
    groups = [collect(Symbol, g) for g in promoters]
    for g in groups, s in g
        addname!(s)
    end
    promoter_set = Set{Symbol}(Iterators.flatten(groups))
    sp = [Species(n, n in promoter_set ? :promoter : :molecule) for n in names]
    sidx = Dict{Symbol,Int}(n => i for (i, n) in enumerate(names))
    N = length(names)
    # --- parameters ----------------------------------------------------------
    pnames = Symbol[]
    literal = Dict{Symbol,Float64}()
    given = Dict{Symbol,Float64}()
    for (k, v) in pairs(params)
        given[Symbol(k)] = Float64(v)
    end
    function register!(x::ParamOrValue, auto::Symbol)
        if x isa Symbol
            x in pnames || push!(pnames, x)
            return x
        else
            auto in pnames || push!(pnames, auto)
            literal[auto] = x
            return auto
        end
    end
    kin = CompiledKinetics[]
    stoich = Vector{Vector{Tuple{Int,Int}}}(undef, length(rxns))
    rstoich = Vector{Vector{Tuple{Int,Int}}}(undef, length(rxns))
    for (j, r) in enumerate(rxns)
        order = reaction_order(r)
        vexp = r.volume == :auto ? float(1 - order) :
               r.volume == :none ? 0.0 :
               r.volume == :proportional ? 1.0 : -1.0
        reac = [(sidx[first(x)], last(x)) for x in r.reactants]
        rstoich[j] = reac
        net = Dict{Int,Int}()
        for (s, c) in r.reactants
            net[sidx[s]] = get(net, sidx[s], 0) - c
        end
        for (s, c) in r.products
            net[sidx[s]] = get(net, sidx[s], 0) + c
        end
        stoich[j] = sort!([(s, c) for (s, c) in net if c != 0])
        tag = _sanitize(r.name)
        k = r.kinetics
        if k isa MassAction
            rate = register!(k.rate, k.rate)
            push!(kin, MassActionC(0, reac, vexp, r.copy_number))
        elseif k isa Hill
            register!(k.rate, k.rate)
            for i in eachindex(k.regulators)
                register!(k.K[i], Symbol("__", tag, "_K", i))
                register!(k.n[i], Symbol("__", tag, "_n", i))
            end
            register!(k.basal, Symbol("__", tag, "_basal"))
            push!(kin, HillC(0, reac, vexp, r.copy_number, Int[], Bool[], Int[], Int[], 0, k.logic == :and, k.concentration))
        elseif k isa Custom
            foreach(s -> register!(s, s), k.params)
            push!(kin, CustomC(k.f, Int[], reac, vexp, r.copy_number))
        else
            throw(ArgumentError("unknown kinetics $(typeof(k))"))
        end
    end
    pidx = Dict{Symbol,Int}(n => i for (i, n) in enumerate(pnames))
    p0 = Vector{Float64}(undef, length(pnames))
    for (i, n) in enumerate(pnames)
        if haskey(literal, n)
            p0[i] = literal[n]
        elseif haskey(given, n)
            p0[i] = given[n]
        else
            throw(ArgumentError("no value given for parameter :$n (pass it via `params`)"))
        end
    end
    # second pass: resolve parameter indices
    for (j, r) in enumerate(rxns)
        k = r.kinetics
        tag = _sanitize(r.name)
        c = kin[j]
        if k isa MassAction
            kin[j] = MassActionC(pidx[k.rate], c.reactants, c.vol_exp, c.copy_number)
        elseif k isa Hill
            regs = [sidx[s] for s in k.regulators]
            act = [m == :activate for m in k.modes]
            Kidx = [pidx[k.K[i] isa Symbol ? k.K[i] : Symbol("__", tag, "_K", i)] for i in eachindex(regs)]
            nidx = [pidx[k.n[i] isa Symbol ? k.n[i] : Symbol("__", tag, "_n", i)] for i in eachindex(regs)]
            bidx = pidx[k.basal isa Symbol ? k.basal : Symbol("__", tag, "_basal")]
            kin[j] = HillC(pidx[k.rate], c.reactants, c.vol_exp, c.copy_number, regs, act, Kidx, nidx, bidx, c.and_logic, c.conc)
        else
            kin[j] = CustomC(k.f, [pidx[s] for s in k.params], c.reactants, c.vol_exp, c.copy_number)
        end
    end
    # --- dependency graph ------------------------------------------------------
    deps = Vector{Set{Int}}(undef, length(rxns))
    for (j, r) in enumerate(rxns)
        d = Set{Int}(first.(rstoich[j]))
        c = kin[j]
        c isa HillC && union!(d, c.regs)
        c isa CustomC && union!(d, 1:N)
        deps[j] = d
    end
    depgraph = Vector{Vector{Int}}(undef, length(rxns))
    for j in eachindex(rxns)
        changed = Set{Int}(first.(stoich[j]))
        depgraph[j] = [k for k in eachindex(rxns) if !isempty(intersect(deps[k], changed))]
    end
    # --- highest order of reaction per species (for adaptive tau) ---------------
    hor = zeros(Int, N)
    mrr = zeros(Int, N)
    for j in eachindex(rxns)
        order = sum(last, rstoich[j]; init=0)
        for (s, c) in rstoich[j]
            hor[s] = max(hor[s], order)
            mrr[s] = max(mrr[s], c)
        end
    end
    pg = [[sidx[s] for s in g] for g in groups]
    isprom = [s.kind == :promoter for s in sp]
    ReactionModel(sp, sidx, rxns, kin, stoich, rstoich, depgraph, pnames, pidx, p0, hor, mrr, pg, isprom)
end

ReactionModel(reactions::Reaction...; kwargs...) = ReactionModel(collect(reactions); kwargs...)

nspecies(m::ReactionModel) = length(m.species)
nreactions(m::ReactionModel) = length(m.reactions)
nparams(m::ReactionModel) = length(m.pnames)
speciesnames(m::ReactionModel) = [s.name for s in m.species]
paramnames(m::ReactionModel) = copy(m.pnames)
speciesindex(m::ReactionModel, s::Symbol) = m.sidx[s]
speciesindex(m::ReactionModel, i::Integer) = Int(i)
paramindex(m::ReactionModel, s::Symbol) = m.pidx[s]

"""
    stoichiometry(model) -> Matrix{Int}

Dense net stoichiometry matrix (species × reactions).
"""
function stoichiometry(m::ReactionModel)
    S = zeros(Int, nspecies(m), nreactions(m))
    for j in 1:nreactions(m), (s, c) in m.stoich[j]
        S[s, j] = c
    end
    S
end

"""
    set_params(model, p=model.p0; kwargs...) -> Vector{Float64}

Return a copy of the parameter vector `p` with the named overrides applied,
e.g. `set_params(model; k_on = 0.5)`.
"""
function set_params(m::ReactionModel, p::AbstractVector{<:Real}=m.p0; kwargs...)
    q = Vector{Float64}(p)
    for (k, v) in kwargs
        haskey(m.pidx, k) || throw(ArgumentError("unknown parameter :$k"))
        q[m.pidx[k]] = Float64(v)
    end
    q
end
set_params(m::ReactionModel, d::AbstractDict) = set_params(m; (Symbol(k) => v for (k, v) in d)...)

"""
    initial_state(model; kwargs...) -> Vector{Int}
    initial_state(model, dict)

Integer count vector with the given species values (unlisted species are 0),
e.g. `initial_state(model; mRNA = 5, G_off = 1)`.
"""
function initial_state(m::ReactionModel; kwargs...)
    x = zeros(Int, nspecies(m))
    for (k, v) in kwargs
        haskey(m.sidx, k) || throw(ArgumentError("unknown species :$k"))
        x[m.sidx[k]] = Int(v)
    end
    x
end
initial_state(m::ReactionModel, d::AbstractDict) = initial_state(m; (Symbol(k) => v for (k, v) in d)...)

function Base.show(io::IO, m::ReactionModel)
    print(io, "ReactionModel(", nspecies(m), " species, ", nreactions(m), " reactions, ", nparams(m), " parameters)")
end
function Base.show(io::IO, ::MIME"text/plain", m::ReactionModel)
    println(io, "ReactionModel with ", nspecies(m), " species, ", nreactions(m), " reactions, ", nparams(m), " parameters")
    println(io, "  species: ", join(("$(s.name)$(s.kind == :promoter ? "*" : "")" for s in m.species), ", "))
    for r in m.reactions
        println(io, "  ", r)
    end
    print(io, "  parameters: ", join(("$(n)=$(round(m.p0[i], sigdigits=4))" for (i, n) in enumerate(m.pnames) if !startswith(String(n), "__")), ", "))
end
