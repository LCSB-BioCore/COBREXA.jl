using Documenter, CobraTools


makedocs(
    sitename="CobraTools.jl",
    authors = "St. Elmo Wilken",
    pages = [
        "Home" => "index.md",
        "Model IO" => "io.md",
        "Model Construction" => "model_construction.md",
        "Optimization Based Analysis" => "basic_analysis.md",
        "Sampling Tools" => "sampling_tools.md",
        "Equilibrator Interface" => "thermo_tools.md",
        "Brenda Interface" => "brenda_tools.md"
    ]
)

deploydocs(
    repo = "github.com/stelmo/CobraTools.jl.git",
)