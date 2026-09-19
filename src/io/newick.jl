"""
    newick(lineage_table; founder=nothing, t_end=nothing) -> String

Newick representation of the division tree(s). Branch lengths are cell
lifetimes; cells still alive get length `t_end - birth_time` when `t_end` is
given. With `founder` only that founder's tree is returned; otherwise all trees
are joined under a virtual root.
"""
function newick(lt::LineageTable; founder=nothing, t_end=nothing)
    cm = children_map(lt)
    function len(id)
        e = lt.end_time[id]
        isnan(e) && (e = t_end === nothing ? lt.birth_time[id] : t_end)
        max(e - lt.birth_time[id], 0.0)
    end
    function node(id)
        ch = get(cm, id, Int[])
        inner = isempty(ch) ? "" : "(" * join((node(c) for c in ch), ",") * ")"
        string(inner, "c", id, ":", round(len(id), digits = 6))
    end
    if founder !== nothing
        return node(Int(founder)) * ";"
    end
    founders = lt.id[lt.parent .== 0]
    length(founders) == 1 && return node(founders[1]) * ";"
    "(" * join((node(f) for f in founders), ",") * ")root;"
end
