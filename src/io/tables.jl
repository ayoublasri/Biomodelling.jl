# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

"""
    write_counts_csv(path, counts, genes, cell_ids)

Write a cells × genes count matrix with a header row (`cell_id`, gene names).
"""
function write_counts_csv(path::AbstractString, counts::AbstractMatrix, genes::AbstractVector, cell_ids::AbstractVector)
    open(path, "w") do io
        println(io, join(vcat(["cell_id"], string.(genes)), ","))
        for i in 1:size(counts, 1)
            println(io, string(cell_ids[i]), ",", join(string.(view(counts, i, :)), ","))
        end
    end
    path
end

"""
    write_metadata_csv(path, snapshot)

Write per-cell metadata (volume, age, generation, clone, copies, perturbed) of a
snapshot returned by [`snapshot`](@ref) or [`sample_cells`](@ref).
"""
function write_metadata_csv(path::AbstractString, sn::NamedTuple)
    open(path, "w") do io
        println(io, "cell_id,volume,age,generation,clone,copies,perturbed")
        for i in eachindex(sn.ids)
            println(io, sn.ids[i], ",", sn.volume[i], ",", sn.age[i], ",", sn.generation[i], ",", sn.clone[i], ",",
                    sn.copies[i], ",", Int(sn.perturbed[i]))
        end
    end
    path
end

"""
    write_lineage_csv(path, lineage_table)
"""
function write_lineage_csv(path::AbstractString, lt::LineageTable)
    open(path, "w") do io
        println(io, "id,parent,clone,generation,birth_time,end_time,fate,V_birth,V_end")
        for i in eachindex(lt.id)
            println(io, lt.id[i], ",", lt.parent[i], ",", lt.clone[i], ",", lt.generation[i], ",", lt.birth_time[i], ",",
                    lt.end_time[i], ",", lt.fate[i], ",", lt.V_birth[i], ",", lt.V_end[i])
        end
    end
    path
end

"""
    write_h5ad(path, counts; obs, var, obsm, uns)

Write an AnnData-compatible HDF5 file. Requires `using HDF5` (package
extension).
"""
function write_h5ad(args...; kwargs...)
    error("write_h5ad requires the HDF5 package: run `using HDF5` before calling it")
end
