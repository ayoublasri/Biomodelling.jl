abstract type Partitioning end

"""
    BinomialPartition(; σ=0.0)

At division the volume fraction inherited by the first daughter is
`f ~ Normal(0.5, σ)` (clipped to [0.05, 0.95]) and every molecule is assigned to
that daughter independently with probability `f` (binomial partitioning).
Promoter species are inherited, not partitioned (see [`ReactionModel`](@ref)).
"""
struct BinomialPartition <: Partitioning
    σ::Float64
end
BinomialPartition(; σ::Real=0.0) = BinomialPartition(Float64(σ))

"""
    BetaBinomialPartition(; σ=0.0, ρ=0.05)

Over-dispersed partitioning: the per-molecule probability is drawn from a Beta
distribution with mean `f` and intra-class correlation `ρ`, then molecules are
partitioned binomially (models clustered or aggregated molecules).
"""
struct BetaBinomialPartition <: Partitioning
    σ::Float64
    ρ::Float64
end
BetaBinomialPartition(; σ::Real=0.0, ρ::Real=0.05) = BetaBinomialPartition(Float64(σ), Float64(ρ))

partition_fraction(pt::Partitioning, rng::AbstractRNG) = clamp(0.5 + pt.σ * randn(rng), 0.05, 0.95)

@inline function partition_count(n::Int, f::Float64, ::BinomialPartition, rng::AbstractRNG)
    n == 0 && return 0
    rand(rng, Binomial(n, f))
end
@inline function partition_count(n::Int, f::Float64, pt::BetaBinomialPartition, rng::AbstractRNG)
    n == 0 && return 0
    α = f * (1 - pt.ρ) / pt.ρ
    β = (1 - f) * (1 - pt.ρ) / pt.ρ
    rand(rng, Binomial(n, rand(rng, Beta(α, β))))
end

"""
    partition!(xa, xb, x, f, model, pt, base_copies, rng)

Split the state `x` of a dividing cell into daughters `xa` and `xb`. Molecules
are partitioned according to `pt`; promoter groups are inherited by both
daughters when unreplicated (total count equals `base_copies`) and split into
equal halves at random when replicated.
"""
function partition!(xa::Vector{Int}, xb::Vector{Int}, x::Vector{Int}, f::Float64, m::ReactionModel,
                    pt::Partitioning, base_copies::Vector{Int}, rng::AbstractRNG)
    @inbounds for s in eachindex(x)
        if m.is_promoter[s]
            xa[s] = 0; xb[s] = 0
        else
            k = partition_count(x[s], f, pt, rng)
            xa[s] = k; xb[s] = x[s] - k
        end
    end
    for (gi, g) in enumerate(m.promoter_groups)
        total = sum(x[s] for s in g)
        if total <= base_copies[gi]
            for s in g
                xa[s] = x[s]; xb[s] = x[s]
            end
        else
            states = Int[]
            for s in g, _ in 1:x[s]
                push!(states, s)
            end
            shuffle!(rng, states)
            h = cld(length(states), 2)
            for (i, s) in enumerate(states)
                i <= h ? (xa[s] += 1) : (xb[s] += 1)
            end
        end
    end
    xa, xb
end

function randomize_promoters!(x::Vector{Int}, m::ReactionModel, base_copies::Vector{Int}, rng::AbstractRNG)
    for (gi, g) in enumerate(m.promoter_groups)
        for s in g
            x[s] = 0
        end
        for _ in 1:base_copies[gi]
            x[rand(rng, g)] += 1
        end
    end
    x
end
