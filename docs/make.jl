using Documenter
using Literate, JSON
using COBREXA

# some settings
dev_docs_folder = "dev"
pages_branch = "gh-pages"
github_repo_slug = "LCSB-BioCore/COBREXA.jl"
delete!(ENV, "GITHUB_REPOSITORY")

# generate notebooks
notebooks_path = joinpath(@__DIR__, "src", "notebooks")
notebooks_basenames = filter(x -> endswith(x, ".jl"), readdir(notebooks_path))
@info "base names:" notebooks_basenames
notebooks = joinpath.(notebooks_path, notebooks_basenames)
notebooks_outdir = joinpath(@__DIR__, "src", "notebooks")



for notebook in notebooks
    Literate.markdown(
        notebook,
        notebooks_outdir;
        repo_root_url = "https://github.com/$github_repo_slug/blob/master",
        nbviewer_root_url = "https://nbviewer.jupyter.org/github/$github_repo_slug/blob/gh-pages/$dev_docs_folder",
        binder_root_url = "https://mybinder.org/v2/gh/$github_repo_slug/$pages_branch?filepath=$dev_docs_folder",
    )
    Literate.notebook(notebook, notebooks_outdir)
end

# generate index.md from .template and the quickstart in README.md
readme = open(f -> read(f, String), joinpath(@__DIR__, "..", "README.md"))
quickstart =
    match(r"<!--quickstart_begin-->\n([^\0]*)<!--quickstart_end-->", readme).captures[1]
acks = match(
    r"<!--acknowledgements_begin-->\n([^\0]*)<!--acknowledgements_end-->",
    readme,
).captures[1]
ack_logos =
    match(r"<!--ack_logos_begin-->\n([^\0]*)<!--ack_logos_end-->", readme).captures[1]
index_md = open(f -> read(f, String), joinpath(@__DIR__, "src", "index.md.template"))
index_md = replace(index_md, "<!--insert_quickstart-->\n" => quickstart)
index_md = replace(index_md, "<!--insert_acknowledgements-->\n" => acks)
index_md = replace(index_md, "<!--insert_ack_logos-->\n" => ack_logos)
open(f -> write(f, index_md), joinpath(@__DIR__, "src", "index.md"), "w")

# copy the contribution guide
cp(
    joinpath(@__DIR__, "..", ".github", "CONTRIBUTING.md"),
    joinpath(@__DIR__, "src", "howToContribute.md"),
    force = true,
)

find_mds(path) =
    joinpath.(
        Ref(path),
        filter(x -> endswith(x, ".md"), readdir(joinpath(@__DIR__, "src", path))),
    )

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
        "User guide" => [
            "Quickstart tutorials" => vcat(
                "Detailed tutorial listing" => "tutorials.md",
                find_mds("tutorials"),
            ),
            "Advanced tutorials" => vcat(
                "Detailed tutorial listing" => "advanced.md",
                find_mds("advanced"),
            ),
            "Examples and notebooks" => vcat(
                "Detailed notebook listing" => "notebooks.md",
                find_mds("notebooks"),
            ),
        ],
        "API reference" => vcat("Contents" => "functions.md", find_mds("functions")),
        "How to contribute" => "howToContribute.md",
    ],
)

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

# clean up notebooks -- we do not need to deploy all the stuff that was
# generated in the process
#
# extra fun: failing programs (such as plotting libraries) may generate core
# dumps that contain the dumped environment strings, which in turn contain
# github auth tokens. These certainly need to be avoided.
notebooks_names = [n[begin:end-3] for n in notebooks_basenames]
ipynb_names = notebooks_names .* ".ipynb"
notebooks_allowed_files = vcat("index.html", ipynb_names)
@info "allowed files:" notebooks_allowed_files
for (root, dirs, files) in walkdir(joinpath(@__DIR__, "build", "notebooks"))
    for f in files
        if !(f in notebooks_allowed_files)
            @info "removing notebook build artifact `$(joinpath(root, f))'"
            rm(joinpath(root, f))
        end
    end
end

# also remove the index template
rm(joinpath(@__DIR__, "build", "index.md.template"))

# Binder actually has 1.6.2 kernel (seen in October 2021), but for whatever
# reason it's called julia-1.1. Don't ask me.
for ipynb in joinpath.(@__DIR__, "build", "notebooks", ipynb_names)
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
