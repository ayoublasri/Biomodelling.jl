module BiomodellingHDF5Ext

using Biomodelling, HDF5

function _write_array!(g, name::AbstractString, A)
    if A isa AbstractMatrix
        ds = write_dataset(g, name, permutedims(Array(A)))
    else
        ds = write_dataset(g, name, Array(A))
    end
    attrs(g[name])["encoding-type"] = "array"
    attrs(g[name])["encoding-version"] = "0.2.0"
end

function _write_column!(g, name::AbstractString, col)
    if eltype(col) <: AbstractString || eltype(col) <: Symbol
        write_dataset(g, name, String.(col))
        attrs(g[name])["encoding-type"] = "string-array"
        attrs(g[name])["encoding-version"] = "0.2.0"
    elseif eltype(col) <: Bool
        write_dataset(g, name, Int8.(col))
        attrs(g[name])["encoding-type"] = "array"
        attrs(g[name])["encoding-version"] = "0.2.0"
    else
        write_dataset(g, name, Array(col))
        attrs(g[name])["encoding-type"] = "array"
        attrs(g[name])["encoding-version"] = "0.2.0"
    end
end

function _write_dataframe!(parent, name::AbstractString, index::AbstractVector, cols)
    g = create_group(parent, name)
    attrs(g)["encoding-type"] = "dataframe"
    attrs(g)["encoding-version"] = "0.2.0"
    attrs(g)["_index"] = "_index"
    names = String[]
    for (k, v) in pairs(cols)
        push!(names, String(k))
        _write_column!(g, String(k), v)
    end
    attrs(g)["column-order"] = names
    write_dataset(g, "_index", String.(index))
    attrs(g["_index"])["encoding-type"] = "string-array"
    attrs(g["_index"])["encoding-version"] = "0.2.0"
    g
end

"""
    write_h5ad(path, counts; genes=nothing, cells=nothing, obs=(;), var=(;), obsm=(;), uns=(;))

Write a cells × genes matrix and metadata as an AnnData-compatible `.h5ad` file
(readable with `anndata.read_h5ad`).
"""
function Biomodelling.write_h5ad(path::AbstractString, counts::AbstractMatrix{<:Real}; genes=nothing, cells=nothing,
                                 obs=NamedTuple(), var=NamedTuple(), obsm=NamedTuple(), uns=NamedTuple())
    n, g = size(counts)
    cellnames = cells === nothing ? ["cell_$i" for i in 1:n] : string.(cells)
    genenames = genes === nothing ? ["gene_$j" for j in 1:g] : string.(genes)
    h5open(path, "w") do f
        attrs(f)["encoding-type"] = "anndata"
        attrs(f)["encoding-version"] = "0.1.0"
        _write_array!(f, "X", Float32.(counts))
        _write_dataframe!(f, "obs", cellnames, obs)
        _write_dataframe!(f, "var", genenames, var)
        gm = create_group(f, "obsm")
        attrs(gm)["encoding-type"] = "dict"; attrs(gm)["encoding-version"] = "0.1.0"
        for (k, v) in pairs(obsm)
            _write_array!(gm, String(k), v)
        end
        gu = create_group(f, "uns")
        attrs(gu)["encoding-type"] = "dict"; attrs(gu)["encoding-version"] = "0.1.0"
        for (k, v) in pairs(uns)
            if v isa AbstractString || v isa Symbol
                write_dataset(gu, String(k), String(v))
                attrs(gu[String(k)])["encoding-type"] = "string"; attrs(gu[String(k)])["encoding-version"] = "0.2.0"
            elseif v isa Number
                write_dataset(gu, String(k), v)
                attrs(gu[String(k)])["encoding-type"] = "numeric-scalar"; attrs(gu[String(k)])["encoding-version"] = "0.2.0"
            else
                _write_array!(gu, String(k), v)
            end
        end
        for name in ("layers", "obsp", "varm", "varp")
            gg = create_group(f, name)
            attrs(gg)["encoding-type"] = "dict"; attrs(gg)["encoding-version"] = "0.1.0"
        end
    end
    path
end

end # module
