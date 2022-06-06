using Documenter
using Literate, JSON
using COBREXA

# some settings
dev_docs_folder = "dev"
pages_branch = "gh-pages"

# This must match the repo slug on github!
github_repo_slug = ENV["CI_PROJECT_NAMESPACE"] * "/" * ENV["CI_PROJECT_NAME"]

# generate examples
examples_path = joinpath(@__DIR__, "src", "examples")
examples_basenames = filter(x -> endswith(x, ".jl"), readdir(examples_path))
@info "base names:" examples_basenames
examples = joinpath.(examples_path, examples_basenames)
examples_outdir = joinpath(@__DIR__, "src", "examples")

for example in examples
    #TODO improve how the nbviewer and binder links are inserted. Direct link to ipynb would be cool
    Literate.markdown(
        example,
        examples_outdir;
        repo_root_url = "https://github.com/$github_repo_slug/blob/master",
        nbviewer_root_url = "https://nbviewer.jupyter.org/github/$github_repo_slug/blob/gh-pages/$dev_docs_folder",
        binder_root_url = "https://mybinder.org/v2/gh/$github_repo_slug/$pages_branch?filepath=$dev_docs_folder",
    )
    Literate.notebook(example, exampless_outdir)
end

# extract shared documentation parts from README.md
readme_md = open(f -> read(f, String), joinpath(@__DIR__, "..", "README.md"))
quickstart =
    match(r"<!--quickstart_begin-->\n([^\0]*)<!--quickstart_end-->", readme_md).captures[1]
acks = match(
    r"<!--acknowledgements_begin-->\n([^\0]*)<!--acknowledgements_end-->",
    readme_md,
).captures[1]
ack_logos =
    match(r"<!--ack_logos_begin-->\n([^\0]*)<!--ack_logos_end-->", readme_md).captures[1]

# insert the shared documentation parts into index and quickstart templates
#TODO use direct filename read/write
index_md = open(f -> read(f, String), joinpath(@__DIR__, "src", "index.md.template"))
index_md = replace(index_md, "<!--insert_acknowledgements-->\n" => acks)
index_md = replace(index_md, "<!--insert_ack_logos-->\n" => ack_logos)
open(f -> write(f, index_md), joinpath(@__DIR__, "src", "index.md"), "w")

quickstart_md =
    open(f -> read(f, String), joinpath(@__DIR__, "src", "quickstart.md.template"))
quickstart_md = replace(quickstart_md, "<!--insert_quickstart-->\n" => quickstart)
open(f -> write(f, quickstart_md), joinpath(@__DIR__, "src", "quickstart.md"), "w")

# copy the contribution guide
cp(
    joinpath(@__DIR__, "..", ".github", "CONTRIBUTING.md"),
    joinpath(@__DIR__, "src", "howToContribute.md"),
    force = true,
)

# a helper for sourcing the documentation files from directories
find_mds(path) =
    joinpath.(
        Ref(path),
        filter(x -> endswith(x, ".md"), readdir(joinpath(@__DIR__, "src", path))),
    )

# Documenter tries to guess the repo slug from git remote URL but that doesn't
# work really well here, this is the only fallback. If this breaks, "Edit on
# GitHub" links will stop working. (See Documenter.jl source in
# src/Utilities/Utilities.jl, in November 2021 it was around line 500) -mk
ENV["TRAVIS_REPO_SLUG"] = github_repo_slug

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
        "Quick start" => "quickstart.md",
        "User guide" => [
            "Examples and notebooks" =>
                vcat("All examples" => "examples.md", find_mds("examples")),
            "Core concepts" =>
                vcat("All tutorials" => "tutorials.md", find_mds("tutorials")),
        ],
        "Function reference" => vcat("Contents" => "functions.md", find_mds("functions")),
        "How to contribute" => "howToContribute.md",
    ],
)

# remove the workaround (this would cause deploydocs() to get confused and try
# to deploy the travis way)
delete!(ENV, "TRAVIS_REPO_SLUG")

# replace the "edit this" links for the generated documentation
function replace_in_doc(filename, replacement)
    contents = open(f -> read(f, String), joinpath(@__DIR__, "build", filename))
    contents = replace(contents, replacement)
    open(f -> write(f, contents), joinpath(@__DIR__, "build", filename), "w")
end

replace_in_doc("index.html", "blob/master/docs/src/index.md" => "")
replace_in_doc(
    joinpath("howToContribute", "index.html"),
    "blob/master/docs/src/howToContribute.md" => "blob/master/.github/CONTRIBUTING.md",
)

# clean up examples -- we do not need to deploy all the stuff that was
# generated in the process
#
# extra fun: failing programs (such as plotting libraries) may generate core
# dumps that contain the dumped environment strings, which in turn contain
# github auth tokens. These certainly need to be avoided.
examples_names = [n[begin:end-3] for n in examples_basenames]
ipynb_names = examples_names .* ".ipynb"
examples_allowed_files = vcat("index.html", ipynb_names)
@info "allowed files:" examples_allowed_files
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "build", "examples"))
    for f in files
        if !(f in examples_allowed_files)
            @info "removing notebook build artifact `$(joinpath(root, f))'"
            rm(joinpath(root, f))
        end
    end
end

# remove the template files
rm(joinpath(@__DIR__, "build", "index.md.template"))
rm(joinpath(@__DIR__, "build", "quickstart.md.template"))

# Binder actually has 1.6.2 kernel (seen in October 2021), but for whatever
# reason it's called julia-1.1. Don't ask me.
for ipynb in joinpath.(@__DIR__, "build", "examples", ipynb_names)
    @info "changing julia version to 1.1 in `$ipynb'"
    js = JSON.parsefile(ipynb)
    js["metadata"]["kernelspec"]["name"] = "julia-1.1"
    js["metadata"]["kernelspec"]["display_name"] = "Julia 1.1.0"
    js["metadata"]["language_info"]["version"] = "1.1.0"
    open(f -> JSON.print(f, js), ipynb, "w")
end

# deploy the result
deploydocs(
    repo = "github.com/$github_repo_slug.git",
    target = "build",
    branch = pages_branch,
    devbranch = "develop",
)
