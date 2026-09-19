# ---------------------------------------------------------------------------
# Luria-Delbrück style fluctuation tests and clonal memory scores
# ---------------------------------------------------------------------------

"""
    fluctuation_test(model, x0, phenotype; n_clones=50, generations=6, settings, p, rng)

Grow `n_clones` clones from single cells for `generations` doublings under free
growth and evaluate `phenotype(snapshot)` (a fraction in [0, 1], e.g. the
fraction of cells above a threshold) in each clone. Returns the per-clone
values, their mean and variance, the binomial (no-memory) variance expected
for the observed clone sizes, and the variance ratio: values well above 1
indicate heritable states.
"""
function fluctuation_test(m::ReactionModel, x0::AbstractVector{<:Integer}, phenotype::Function;
                          n_clones::Int=50, generations::Int=6, settings::PopulationSettings=PopulationSettings(),
                          p::AbstractVector{<:Real}=m.p0, rng::AbstractRNG=Random.default_rng(), perturbation=nothing)
    T = generations * doubling_time(settings.growth)
    st = PopulationSettings(; dt = settings.dt, kernel = settings.kernel, growth = settings.growth,
                            size_control = settings.size_control, partitioning = settings.partitioning,
                            replication = settings.replication, control = FreeGrowth(max_cells = 2^(generations + 3)),
                            background_death = settings.background_death, record_every = 10^6,
                            track_lineage = false, randomize_promoters = settings.randomize_promoters,
                            threads = settings.threads)
    vals = Float64[]; sizes = Int[]
    for _ in 1:n_clones
        res = simulate_population(m, x0, 1, (0.0, T); settings = st, p = p, rng = Xoshiro(rand(rng, UInt64)), perturbation = perturbation)
        sn = final_snapshot(res)
        push!(vals, Float64(phenotype(sn))); push!(sizes, size(sn.counts, 1))
    end
    μ = mean(vals)
    binom = mean(μ * (1 - μ) ./ max.(sizes, 1))
    (values = vals, sizes = sizes, mean = μ, var = var(vals), binomial_var = binom, ratio = var(vals) / max(binom, eps()))
end

"""
    clonal_variance_scores(result; record=length(result.t), n_perm=200, min_cells=3,
                           normalize=:concentration, rng) -> (score, pvalue, clone_means, clones)
    clonal_variance_scores(Y, clones; n_perm=200, min_cells=3, rng)

MemorySeq-style heritability score per gene: variance of clone-mean expression
divided by its expectation under random clone labels (permutation null).
Scores well above 1 flag heritable ("memory") genes. The second form works on
any cells × genes matrix `Y` (for example observed counts) with clone labels.
"""
function clonal_variance_scores(r::PopulationResult; record::Int=length(r.t), n_perm::Int=200, min_cells::Int=3,
                                normalize::Symbol=:concentration, rng::AbstractRNG=Random.default_rng())
    Y = normalize == :concentration ? r.counts[record] ./ r.volume[record] : float.(r.counts[record])
    clonal_variance_scores(Y, r.clone[record]; n_perm, min_cells, rng)
end

function clonal_variance_scores(Y::AbstractMatrix{<:Real}, cl::AbstractVector; n_perm::Int=200, min_cells::Int=3,
                                rng::AbstractRNG=Random.default_rng())
    clones = [c for (c, n) in countmap(cl) if n >= min_cells]
    length(clones) < 3 && throw(ArgumentError("need at least 3 clones with ≥ $min_cells cells"))
    idx = [findall(==(c), cl) for c in clones]
    G = size(Y, 2)
    function clone_var(labels_idx)
        M = zeros(G)
        for g in 1:G
            means = [mean(view(Y, ii, g)) for ii in labels_idx]
            M[g] = var(means)
        end
        M
    end
    obs = clone_var(idx)
    sizes = length.(idx)
    null = zeros(n_perm, G)
    allrows = reduce(vcat, idx)
    for k in 1:n_perm
        perm = shuffle(rng, allrows)
        pidx = Vector{Vector{Int}}(undef, length(idx))
        pos = 1
        for (q, n) in enumerate(sizes)
            pidx[q] = perm[pos:pos + n - 1]; pos += n
        end
        null[k, :] = clone_var(pidx)
    end
    nullmean = vec(mean(null; dims = 1))
    score = obs ./ max.(nullmean, eps())
    pval = [(1 + count(>=(obs[g]), view(null, :, g))) / (n_perm + 1) for g in 1:G]
    cm = zeros(length(clones), G)
    for (q, ii) in enumerate(idx), g in 1:G
        cm[q, g] = mean(view(Y, ii, g))
    end
    (score = score, pvalue = pval, clone_means = cm, clones = clones)
end

"""
    covariance_eigenspectrum(Y; log_transform=true) -> Vector{Float64}

Descending eigenvalues of the gene-gene covariance matrix of (optionally
log1p-transformed) expression, which shares its non-zero spectrum with the
cell-cell covariance matrix.
"""
function covariance_eigenspectrum(Y::AbstractMatrix{<:Real}; log_transform::Bool=true)
    Z = log_transform ? log1p.(float.(Y)) : float.(Y)
    Z = Z .- mean(Z; dims = 1)
    C = Symmetric(Z' * Z ./ max(size(Z, 1) - 1, 1))
    λ = eigvals(C)
    sort!(λ; rev = true)
    λ[λ .> 1e-12]
end

"""
    powerlaw_tail_exponent(λ; k=10) -> Float64

Slope of `log λ_i` against `log i` over the top `k` eigenvalues (rank-size
exponent); a heavy, power-law-like tail is the signature associated with
memory genes (Ghosh, Chakrabarti and Raju 2025).
"""
function powerlaw_tail_exponent(λ::AbstractVector{<:Real}; k::Int=10)
    k = min(k, length(λ))
    k < 3 && return NaN
    x = log.(1:k); y = log.(λ[1:k])
    xm = mean(x); ym = mean(y)
    sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
end
