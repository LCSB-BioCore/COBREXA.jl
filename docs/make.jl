using Documenter, COBREXA

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
        "Background" => "background.md",
        "Tutorials" => "tutorials.md",
        "Functions" => "functions.md",
        "How to contribute" => "howToContribute.md",
    ],
)
