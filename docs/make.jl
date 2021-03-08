using Documenter, CobraTools

makedocs(
    sitename="CobraTools.jl",
    authors = "St. Elmo Wilken",
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Model Structure" => "model_structure.md",
        "Model IO" => "io.md",
        "Model Construction" => "model_construction.md",
        "Optimization Based Analysis Tools" => "basic_analysis.md",
        "Sampling Tools" => "sampling_tools.md",
        "External Tools" => "external_tools.md",
    ]
)

deploydocs(
    repo = "github.com/stelmo/CobraTools.jl.git",
)