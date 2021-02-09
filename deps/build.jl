using JLD
using Measurements

include("buildfuncs.jl")

@info "Building ΔG files"
for f in readdir(joinpath("..", "data"))
    if contains(f, "kegg_reactions_CC_")
        ph = split(f, "_")[end][1:5]
        mkGibbsDB(joinpath("..", "data", f), joinpath("..", "data", "gibbs_$ph.jld"))
    end
end