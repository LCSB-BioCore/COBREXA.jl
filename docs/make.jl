using Documenter, COBREXA
using Literate

ENV["TRAVIS_REPO_SLUG"] = "LCSB-BioCore/COBREXA.jl"

# generate notebooks
EXAMPLE = joinpath(@__DIR__, "src/notebooks", "example.jl")
OUTPUT = joinpath(@__DIR__, "src/tutorials")

folder = "stable"

## only temporary - will be removed once public
branch = "gh-pages"

Literate.markdown(
    EXAMPLE,
    OUTPUT;
    repo_root_url = "https://github.com/$(ENV["TRAVIS_REPO_SLUG"])/blob/master",
    nbviewer_root_url = "https://nbviewer.jupyter.org/github/$(ENV["TRAVIS_REPO_SLUG"])/blob/gh-pages/$(folder)",
    binder_root_url = "https://mybinder.org/v2/gh/$(ENV["TRAVIS_REPO_SLUG"])/$(branch)?filepath=$(folder)",
)
Literate.notebook(EXAMPLE, OUTPUT)


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
        "Tutorials" => "tutorials.md",
        "Function reference" => "functions.md",
        "How to contribute" => "howToContribute.md",
    ],
)


# deploydocs(
#     repo = "github.com/$(ENV["TRAVIS_REPO_SLUG"])",
#     push_preview=true,
#     deploy_config = deployconfig,
# )
