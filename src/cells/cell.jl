"""
    Cell

State of one cell in a population simulation: identity and ancestry (`id`,
`parent`, `clone`, `generation`, `birth_time`), size (`V`, `V_birth`,
`div_target`), gene copy state (`copies`, `replicated`), growth rate, molecule
counts `x`, the cell-specific parameter vector `p` (after drug modulation), its
own random number generator, and flags `perturbed` and `alive`.
"""
mutable struct Cell
    id::Int
    parent::Int
    clone::Int
    generation::Int
    birth_time::Float64
    V::Float64
    V_birth::Float64
    div_target::Float64
    copies::Int
    replicated::Bool
    growth_rate::Float64
    x::Vector{Int}
    p::Vector{Float64}
    rng::Xoshiro
    perturbed::Bool
    alive::Bool
end

age(c::Cell, t::Real) = t - c.birth_time
Base.show(io::IO, c::Cell) = print(io, "Cell(id=", c.id, ", gen=", c.generation, ", V=", round(c.V, digits=3), ", x=", c.x, ")")
