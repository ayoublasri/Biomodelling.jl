# PIDC network inference (NetworkInference.jl) on the Figure 5b datasets. Run in an environment that
# contains NetworkInference (it pins older dependencies, so it lives in its own environment).
using NetworkInference, DelimitedFiles
const OUT = joinpath(@__DIR__, "..", "output")
for f in filter(x -> startswith(x, "fig5b_data_") && endswith(x, ".csv"), readdir(OUT))
    name = f[length("fig5b_data_")+1:end-4]
    raw = readdlm(joinpath(OUT, f), ',', Float64; skipstart = 1)
    G = size(raw, 2)
    tmp = joinpath(OUT, "tmp_pidc_$name.txt")
    open(tmp, "w") do io
        for j in 1:G
            println(io, "gene_$j\t", join(raw[:, j], "\t"))
        end
    end
    nodes = get_nodes(tmp)
    net = InferredNetwork(PIDCNetworkInference(), nodes)
    M = zeros(G, G)
    for e in net.edges
        i = parse(Int, replace(e.nodes[1].label, "gene_" => "")); j = parse(Int, replace(e.nodes[2].label, "gene_" => ""))
        M[i, j] = e.weight; M[j, i] = e.weight
    end
    open(joinpath(OUT, "fig5b_scores_pidc_$name.csv"), "w") do io
        println(io, join(["gene_$j" for j in 1:G], ","))
        writedlm(io, M, ',')
    end
    rm(tmp)
    println("PIDC done: ", name)
end
