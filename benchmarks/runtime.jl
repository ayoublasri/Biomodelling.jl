# Runtime scaling of simulate_population with the number of cells, genes and kernel.
# Writes benchmarks/runtime.csv
using Biomodelling, Random, Printf
using Random: Xoshiro

function bench(ngenes, ncells, kernel; T = 20.0, dt = 0.1)
    m, _ = random_grn(ngenes; n_activations = ngenes, n_inhibitions = ngenes ÷ 2, telegraph = true, rng = Xoshiro(1))
    x0 = zeros(Int, nspecies(m)); for g in m.promoter_groups; x0[g[1]] = 1; end
    st = PopulationSettings(dt = dt, kernel = kernel, growth = ExponentialGrowth(log(2) / 20), threads = true, track_lineage = false, record_every = 10)
    simulate_population(m, x0, ncells, (0.0, 2.0); settings = st, rng = Xoshiro(2))   # warm-up / compile
    t = @elapsed simulate_population(m, x0, ncells, (0.0, T); settings = st, rng = Xoshiro(2))
    t
end

open(joinpath(@__DIR__, "runtime.csv"), "w") do io
    println(io, "genes,cells,kernel,threads,seconds,cell_updates_per_second")
    for (ng, nc) in ((10, 100), (10, 1000), (10, 10000), (50, 1000), (100, 1000), (200, 1000)),
        (kname, k) in (("DirectSSA", DirectSSA()), ("HybridSSATau", HybridSSATau(0.1)), ("AdaptiveTauLeap", AdaptiveTauLeap()))
        t = bench(ng, nc, k)
        ups = nc * 200 / t
        @printf(io, "%d,%d,%s,%d,%.3f,%.0f\n", ng, nc, kname, Threads.nthreads(), t, ups)
        @printf("%4d genes %6d cells %-16s %7.2fs  %9.0f cell-steps/s\n", ng, nc, kname, t, ups)
    end
end
