# Runs every simulation script of the paper. Outputs go to paper/output/.
const SCRIPTS = ["fig2_engine.jl", "fig3_growth.jl", "fig4_persisters.jl", "fig5_benchmarks.jl", "fig6_inference.jl"]
for s in SCRIPTS
    println("\n==== ", s, " ====")
    t = @elapsed include(joinpath(@__DIR__, "scripts", s))
    println("---- ", s, " done in ", round(t; digits = 1), " s")
end
