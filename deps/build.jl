using JLD
using Measurements
using GZip
using HTTP

include("buildfuncs.jl")

# gibbs_reaction_url(ph::String) = "http://equilibrator.weizmann.ac.il/static/downloads/kegg_reactions_CC_ph$(ph).csv.gz"
# for ph in string.(5:0.5:9) 
#     url = gibbs_reaction_url(ph)
#     r = HTTP.request("GET", gibbs_reaction_url(ph))
#     if r.status == 200
#         in_file_path = joinpath("..", "data", "equilibrator_data.csv.gz")
#         open(in_file_path, "w") do io
#             write(io, r.body)
#         end
#         mkGibbsDB(in_file_path, joinpath("..", "data", "gibbs_$ph.jld"))
#         rm(in_file_path) # clean up    
#     else
#         @warn "Could not download Gibbs data from Equilibrator for pH = $ph"
#     end
# end

@info "Fetching name space data"
metanetx_url = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv"
