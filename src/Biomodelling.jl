module Biomodelling

using Statistics
using StatsBase
using PoissonRandom
using LsqFit
using Random
using DelimitedFiles
using StringDistances
using LinearAlgebra
using SparseArrays

include("src/comp_g.jl")
include("src/comp_hist_1d.jl")
include("src/comp_tau.jl")
include("src/comp_non_tau.jl")
include("src/compute_L.jl")
include("src/non_negative_Poisson_tauleap.jl")
include("src/HO_reaction.jl")
include("src/Donne.jl")
include("src/rates.jl")
include("src/update_rates.jl")
include("src/exponential_growth.jl")
include("src/tauleap.jl")
include("src/tauleapswitch.jl")
include("src/tauleapswitch2.jl")
include("src/ssa.jl")
include("src/ssa2.jl")
include("src/ssa3.jl")
include("src/ssa_switch.jl")
include("src/ssa_switch2.jl")
include("src/ssa_switch3.jl")
include("src/ssa_switch4.jl")
include("src/adaptive_tauleap.jl")
include("src/random_network.jl")

export comp_g
export comp_hist_1d
export comp_tau
export HO_reaction
export Donne
export rates
export update_rates
export exponential_growth
export tauleap
export tauleapswitch
export tauleapswitch2
export ssa
export ssa2
export ssa3
export ssa_switch
export ssa_switch2
export ssa_switch3
export ssa_switch4
export adaptive_tauleap
export comp_non_tau
export compute_L
export non_negative_Poisson_tauleap
export random_network

end