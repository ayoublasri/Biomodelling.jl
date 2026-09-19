"""
    AbstractKernel

Supertype of stochastic simulation kernels: [`DirectSSA`](@ref),
[`TauLeap`](@ref), [`HybridSSATau`](@ref), [`AdaptiveTauLeap`](@ref).
"""
abstract type AbstractKernel end

"""
    Workspace(model)

Preallocated buffers used by the kernels (propensities, leap counts, scratch
state). One workspace per thread is enough.
"""
mutable struct Workspace
    a::Vector{Float64}
    xtmp::Vector{Int}
    counts::Vector{Int}
    mu::Vector{Float64}
    sig::Vector{Float64}
    critical::Vector{Bool}
end
Workspace(m::ReactionModel) = Workspace(zeros(nreactions(m)), zeros(Int, nspecies(m)),
                                        zeros(Int, nreactions(m)), zeros(nspecies(m)),
                                        zeros(nspecies(m)), falses(nreactions(m)))
