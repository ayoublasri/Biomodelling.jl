# ---------------------------------------------------------------------------
# Observation models
# ---------------------------------------------------------------------------

"""
    SeqProtocol(; capture=0.1, capture_cv=0.3, depth=nothing, dropout=0.0,
                 batch_cv=0.0, n_batches=1)

Single-cell RNA-seq observation model: each cell captures every molecule
independently with a cell-specific efficiency drawn from a Beta distribution of
mean `capture` and coefficient of variation `capture_cv`; optionally the reads of
each cell are resampled to a Poisson-distributed total of mean `depth`
(multinomial down-sampling); an extra Bernoulli `dropout` zeroes gene counts;
`batch_cv` scales the capture efficiency per batch.
"""
Base.@kwdef struct SeqProtocol
    capture::Float64 = 0.1
    capture_cv::Float64 = 0.3
    depth::Union{Nothing,Float64} = nothing
    dropout::Float64 = 0.0
    batch_cv::Float64 = 0.0
    n_batches::Int = 1
end

function _beta_from_mean_cv(m::Float64, cv::Float64)
    v = (cv * m)^2
    v >= m * (1 - m) && throw(ArgumentError("capture_cv too large for mean capture $m"))
    common = m * (1 - m) / v - 1
    Beta(m * common, (1 - m) * common)
end

"""
    sequence(counts, protocol=SeqProtocol(); rng, batches=nothing) -> (Y, capture, batch)

Apply the observation model to a cells × genes matrix of true molecule counts.
"""
function sequence(counts::AbstractMatrix{<:Integer}, proto::SeqProtocol=SeqProtocol();
                  rng::AbstractRNG=Random.default_rng(), batches=nothing)
    n, g = size(counts)
    β = proto.capture_cv > 0 ? rand(rng, _beta_from_mean_cv(proto.capture, proto.capture_cv), n) : fill(proto.capture, n)
    b = batches === nothing ? rand(rng, 1:proto.n_batches, n) : Vector{Int}(batches)
    bf = proto.batch_cv > 0 ? exp.(proto.batch_cv .* randn(rng, proto.n_batches) .- proto.batch_cv^2 / 2) : ones(proto.n_batches)
    Y = zeros(Int, n, g)
    for i in 1:n
        e = clamp(β[i] * bf[b[i]], 0.0, 1.0)
        for j in 1:g
            c = counts[i, j]
            c > 0 && (Y[i, j] = rand(rng, Binomial(c, e)))
        end
        if proto.depth !== nothing
            tot = sum(view(Y, i, :))
            if tot > 0
                target = rand(rng, Poisson(proto.depth))
                if target < tot
                    Y[i, :] = rand(rng, Multinomial(target, view(Y, i, :) ./ tot))
                end
            end
        end
        if proto.dropout > 0
            for j in 1:g
                rand(rng) < proto.dropout && (Y[i, j] = 0)
            end
        end
    end
    (Y = Y, capture = β, batch = b)
end

"""
    smfish(counts; efficiency=0.95, rng) -> Matrix{Int}

Single-molecule FISH style observation: independent detection of each molecule
with probability `efficiency`, no cell-specific capture noise.
"""
function smfish(counts::AbstractMatrix{<:Integer}; efficiency::Real=0.95, rng::AbstractRNG=Random.default_rng())
    Y = similar(Matrix{Int}(undef, size(counts)))
    for i in eachindex(counts)
        c = counts[i]
        Y[i] = c > 0 ? rand(rng, Binomial(c, Float64(efficiency))) : 0
    end
    Y
end

"""
    sample_cells(result, record=length(result.t); n=nothing, rng) -> NamedTuple

Random subset of `n` cells from a recorded snapshot (counts and metadata).
"""
function sample_cells(r::PopulationResult, record::Int=length(r.t); n=nothing, rng::AbstractRNG=Random.default_rng())
    sn = snapshot(r, record)
    N = size(sn.counts, 1)
    idx = n === nothing || n >= N ? collect(1:N) : sort!(sample(rng, 1:N, Int(n); replace = false))
    (t = sn.t, counts = sn.counts[idx, :], volume = sn.volume[idx], ids = sn.ids[idx], age = sn.age[idx],
     generation = sn.generation[idx], clone = sn.clone[idx], copies = sn.copies[idx], perturbed = sn.perturbed[idx],
     species = sn.species)
end

library_size(Y::AbstractMatrix) = vec(sum(Y; dims = 2))
