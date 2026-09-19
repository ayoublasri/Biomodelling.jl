using Documenter, Biomodelling

makedocs(
    sitename = "Biomodelling.jl",
    modules = [Biomodelling],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Tutorials" => ["Single cells" => "tutorial_single_cell.md",
                        "Growing and dividing populations" => "tutorial_population.md",
                        "Drug treatment and persisters" => "tutorial_drug.md",
                        "Synthetic single-cell data" => "tutorial_observation.md",
                        "Inference" => "tutorial_inference.md"],
        "Theory notes" => "theory.md",
        "Migration from v1" => "migration.md",
        "API" => "api.md",
    ],
    warnonly = true,
)

deploydocs(repo = "github.com/ayoublasri/Biomodelling.jl.git", devbranch = "master")
