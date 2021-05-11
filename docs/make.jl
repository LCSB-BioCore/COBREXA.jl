using Documenter, COBREXA
using Literate

ENV["TRAVIS_REPO_SLUG"] = "LCSB-BioCore/COBREXA.jl"

# generate notebooks
notebooks_path = joinpath(@__DIR__, "src", "notebooks-src")
notebooks =
    joinpath.(notebooks_path, filter(x -> endswith(x, ".jl"), readdir(notebooks_path)))
notebooks_outdir = joinpath(@__DIR__, "src", "notebooks")

folder = "stable"

## only temporary - will be removed once public
branch = "gh-pages"

for notebook in notebooks
    Literate.markdown(
        notebook,
        notebooks_outdir;
        repo_root_url = "https://github.com/$(ENV["TRAVIS_REPO_SLUG"])/blob/master",
        nbviewer_root_url = "https://nbviewer.jupyter.org/github/$(ENV["TRAVIS_REPO_SLUG"])/blob/gh-pages/$(folder)",
        binder_root_url = "https://mybinder.org/v2/gh/$(ENV["TRAVIS_REPO_SLUG"])/$(branch)?filepath=$(folder)",
    )
    Literate.notebook(notebook, notebooks_outdir)
end

# generate index.md from .template and the quickstart in README.md
quickstart = match(
    r"<!--quickstart_begin-->\n([^\0]*)<!--quickstart_end-->",
    open(f -> read(f, String), joinpath(@__DIR__, "..", "README.md")),
).captures[1]
index_md = replace(
    open(f -> read(f, String), joinpath(@__DIR__, "src", "index-template.md")),
    "<!--insert_quickstart-->\n" => quickstart,
)
open(f -> write(f, index_md), joinpath(@__DIR__, "src", "index.md"), "w")

# build the docs
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
        "Examples and notebooks" => "notebooks.md",
        "Function reference" => "functions.md",
        "How to contribute" => "howToContribute.md",
    ],
)

deploydocs(
    repo = "github.com/$(ENV["TRAVIS_REPO_SLUG"]).git",
    target = "build",
    branch = "gh-pages",
    push_preview = true
)
