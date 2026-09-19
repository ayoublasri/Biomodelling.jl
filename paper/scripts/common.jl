using Biomodelling, Random, Statistics, Distributions, DelimitedFiles, Printf, StatsBase
using Random: Xoshiro
const OUT = joinpath(@__DIR__, "..", "output")
mkpath(OUT)
save_csv(name, header, data) = open(joinpath(OUT, name), "w") do io
    println(io, join(header, ","))
    writedlm(io, data, ',')
end
save_kv(name, pairs) = open(joinpath(OUT, name), "w") do io
    println(io, "key,value")
    for (k, v) in pairs
        println(io, k, ",", v)
    end
end
ks_discrete(sample, pmf) = maximum(abs.([count(<=(k), sample) / length(sample) for k in 0:length(pmf)-1] .- cumsum(pmf)))
empirical_pmf(sample, nmax) = [count(==(n), sample) / length(sample) for n in 0:nmax]
