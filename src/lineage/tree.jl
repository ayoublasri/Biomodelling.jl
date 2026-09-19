"""
    LineageTable

Record of every cell that existed in a population simulation: `id`, `parent`
(0 for founders), `clone` (founder id), `generation`, `birth_time`, `end_time`,
`fate` (`:alive`, `:divided`, `:died`, `:removed`), volume and state at birth
and at the end of the cell's life (`V_birth`, `V_end`, `x_birth`, `x_end`), and
the `perturbed` flag. Row `i` holds the cell with `id == i`.
"""
mutable struct LineageTable
    id::Vector{Int}
    parent::Vector{Int}
    clone::Vector{Int}
    generation::Vector{Int}
    birth_time::Vector{Float64}
    end_time::Vector{Float64}
    fate::Vector{Symbol}
    V_birth::Vector{Float64}
    V_end::Vector{Float64}
    x_birth::Vector{Vector{Int}}
    x_end::Vector{Vector{Int}}
    perturbed::Vector{Bool}
end
LineageTable() = LineageTable(Int[], Int[], Int[], Int[], Float64[], Float64[], Symbol[], Float64[], Float64[],
                              Vector{Int}[], Vector{Int}[], Bool[])

Base.length(lt::LineageTable) = length(lt.id)
Base.show(io::IO, lt::LineageTable) = print(io, "LineageTable(", length(lt), " cells, ",
    count(==(:divided), lt.fate), " divisions, ", count(==(:died), lt.fate), " deaths)")

function record_birth!(lt::LineageTable, c::Cell)
    push!(lt.id, c.id); push!(lt.parent, c.parent); push!(lt.clone, c.clone); push!(lt.generation, c.generation)
    push!(lt.birth_time, c.birth_time); push!(lt.end_time, NaN); push!(lt.fate, :alive)
    push!(lt.V_birth, c.V); push!(lt.V_end, NaN); push!(lt.x_birth, copy(c.x)); push!(lt.x_end, Int[])
    push!(lt.perturbed, c.perturbed)
    lt
end

function record_end!(lt::LineageTable, id::Int, t::Float64, fate::Symbol, x::Vector{Int}, V::Float64)
    lt.end_time[id] = t
    lt.fate[id] = fate
    lt.V_end[id] = V
    lt.x_end[id] = copy(x)
    lt
end

"""
    children(lt, id) -> Vector{Int}
"""
children(lt::LineageTable, id::Integer) = lt.id[lt.parent .== id]

function children_map(lt::LineageTable)
    d = Dict{Int,Vector{Int}}()
    for (i, p) in enumerate(lt.parent)
        p == 0 && continue
        push!(get!(d, p, Int[]), lt.id[i])
    end
    d
end

"""
    sister_pairs(lt) -> Vector{Tuple{Int,Int}}

Pairs of daughter ids sharing the same mother.
"""
function sister_pairs(lt::LineageTable)
    out = Tuple{Int,Int}[]
    for (p, ch) in children_map(lt)
        length(ch) == 2 && push!(out, (ch[1], ch[2]))
    end
    sort!(out)
end

"""
    mother_daughter_pairs(lt) -> Vector{Tuple{Int,Int}}
"""
function mother_daughter_pairs(lt::LineageTable)
    out = Tuple{Int,Int}[]
    for (i, p) in enumerate(lt.parent)
        p == 0 && continue
        push!(out, (p, lt.id[i]))
    end
    out
end

"""
    lineage_of(lt, id) -> Vector{Int}

Ancestors of a cell, starting with the cell itself and ending with its founder.
"""
function lineage_of(lt::LineageTable, id::Integer)
    out = Int[id]
    p = lt.parent[id]
    while p != 0
        push!(out, p)
        p = lt.parent[p]
    end
    out
end

"""
    cousin_pairs(lt) -> Vector{Tuple{Int,Int}}

Pairs of cells whose mothers are sisters.
"""
function cousin_pairs(lt::LineageTable)
    cm = children_map(lt)
    out = Tuple{Int,Int}[]
    for (a, b) in sister_pairs(lt)
        ca = get(cm, a, Int[]); cb = get(cm, b, Int[])
        for i in ca, j in cb
            push!(out, (i, j))
        end
    end
    out
end

"""
    kin_pairs(lt, g) -> Vector{Tuple{Int,Int}}

Pairs of cells whose most recent common ancestor is `g` generations back:
`g = 1` sisters, `g = 2` first cousins, `g = 3` second cousins, and so on.
Each pair is listed once.
"""
function kin_pairs(lt::LineageTable, g::Integer)
    g >= 1 || throw(ArgumentError("g must be at least 1"))
    cm = children_map(lt)
    descendants(id, depth) = depth == 0 ? Int[id] : reduce(vcat, (descendants(c, depth - 1) for c in get(cm, id, Int[])); init = Int[])
    out = Tuple{Int,Int}[]
    for (a, b) in sister_pairs(lt)
        da = descendants(a, g - 1); db = descendants(b, g - 1)
        for i in da, j in db
            push!(out, (i, j))
        end
    end
    out
end
