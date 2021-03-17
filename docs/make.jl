using Documenter, COBREXA

# download core model
download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

makedocs(
    modules = [COBREXA],
    clean = false,
    sitename = "COBREXA.jl",
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        assets = ["assets/favicon.ico"],
        highlights = ["yaml"],
    ),
    authors = "The developers of COBREXA.jl",
    linkcheck = !("skiplinks" in ARGS),
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md",
        "How to contribute" => "howToContribute.md",
        "Model Structure" => "model_structure.md",
        "Model IO" => "io.md",
        "Model Construction" => "model_construction.md",
        "Optimization Based Analysis Tools" => "basic_analysis.md",
        "Sampling Tools" => "sampling_tools.md",
    ],
)

# delete core model
rm("e_coli_core.json")
