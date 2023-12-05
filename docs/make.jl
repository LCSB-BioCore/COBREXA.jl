using Documenter
using Literate, JSON
using COBREXA

# build the examples
examples_path = joinpath(@__DIR__, "src", "examples")
examples_basenames = sort(filter(x -> endswith(x, ".jl"), readdir(examples_path)))
@info "base names:" examples_basenames
examples = joinpath.(examples_path, examples_basenames)
examples_outdir = joinpath(@__DIR__, "src", "examples")

for example in examples
    Literate.markdown(
        example,
        examples_outdir;
        repo_root_url = "https://github.com/LCSB-BioCore/COBREXA.jl/blob/master",
    )
    Literate.notebook(example, examples_outdir)
end

# a helper for sourcing the documentation files from directories
find_mds(path) =
    joinpath.(
        Ref(path),
        filter(x -> endswith(x, ".md"), readdir(joinpath(@__DIR__, "src", path))),
    )

# build the docs
makedocs(
    modules = [COBREXA],
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
        "Examples" => [
            "Contents" => "examples.md"
            find_mds("examples")
        ],
        "Parallel, distributed and HPC processing" => [
            "Contents" => "distributed.md"
            find_mds("distributed")
        ],
        "Core concepts guide" => [
            "Contents" => "concepts.md"
            find_mds("concepts")
        ],
        "Reference" => [
            "Contents" => "reference.md"
            find_mds("reference")
        ],
    ],
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

# deploy the result
deploydocs(
    repo = "github.com/LCSB-BioCore/COBREXA.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
)
