module tauleaping_adaptive

using Distributions
using PoissonRandom

include("comp_g.jl")
include("comp_hist_1d.jl")
include("comp_tau.jl")
include("comp_non_tau.jl")
include("compute_L.jl")
include("non_negative_Poisson_tauleap.jl")
include("HO_reaction.jl")
include("Donne.jl")
include("rates.jl")
include("update_rates.jl")

include("tauleap.jl")
include("tauleapswitch.jl")
include("ssa.jl")
include("ssa_switch.jl")
include("ssa_switch2.jl")
include("ssa_switch3.jl")
include("adaptive_tauleap.jl")

export comp_g
export comp_hist_1d
export comp_tau
export HO_reaction
export Donne
export rates
export update_rates
export tauleap
export tauleapswitch
export ssa
export ssa_switch
export ssa_switch2
export ssa_switch3
export adaptive_tauleap
export comp_non_tau
export compute_L
export non_negative_Poisson_tauleap

end
