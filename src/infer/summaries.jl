# ---------------------------------------------------------------------------
# Summary statistics and distances
# ---------------------------------------------------------------------------

"""
    moment_summaries(Y) -> Vector{Float64}

Per-gene summary statistics of a cells × genes matrix (or a vector for one gene):
mean, coefficient of variation squared, Fano factor and zero fraction,
concatenated gene-wise.
"""
function moment_summaries(Y::AbstractMatrix{<:Real})
    g = size(Y, 2)
    out = zeros(4g)
    for j in 1:g
        y = view(Y, :, j)
        m = mean(y); v = var(y)
        out[j] = m
        out[g + j] = m > 0 ? v / m^2 : 0.0
        out[2g + j] = m > 0 ? v / m : 0.0
        out[3g + j] = count(==(0), y) / length(y)
    end
    out
end
moment_summaries(y::AbstractVector{<:Real}) = moment_summaries(reshape(y, :, 1))

"""
    summary_distance(s, s_obs; scale=nothing) -> Float64

Root-mean-square relative distance between summary vectors, scaled by `scale`
(default: `max(|s_obs|, 1e-8)` element-wise).
"""
function summary_distance(s::AbstractVector{<:Real}, s_obs::AbstractVector{<:Real}; scale=nothing)
    sc = scale === nothing ? max.(abs.(s_obs), 1e-8) : scale
    sqrt(mean(((s .- s_obs) ./ sc) .^ 2))
end

"""
    wasserstein1(x, y) -> Float64

Wasserstein-1 distance between two empirical one-dimensional distributions.
"""
function wasserstein1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    q = range(0, 1; length = 201)
    mean(abs.(quantile(x, q) .- quantile(y, q)))
end
