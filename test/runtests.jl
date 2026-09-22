using Test
using Biomodelling
using Random, Statistics, LinearAlgebra, Distributions
using Random: Xoshiro

const BM = Biomodelling

# Kolmogorov-Smirnov statistic against a discrete reference pmf
function ks_discrete(sample::AbstractVector{<:Integer}, pmf::AbstractVector{<:Real})
    n = length(sample)
    nmax = length(pmf) - 1
    emp = [count(<=(k), sample) / n for k in 0:nmax]
    cdf = cumsum(pmf)
    maximum(abs.(emp .- cdf))
end

@testset "Biomodelling.jl" begin
    include("test_model.jl")
    include("test_kernels.jl")
    include("test_population.jl")
    include("test_exact_population.jl")
    include("test_perturbation.jl")
    include("test_optimize.jl")
    include("test_lineage.jl")
    include("test_observation.jl")
    include("test_inference.jl")
    include("test_generators_io.jl")
    include("test_compat.jl")
    include("test_aqua.jl")
end
