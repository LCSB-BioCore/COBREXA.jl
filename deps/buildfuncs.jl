# """
# mkGibbsDB(equilibrator_file_loc::String)

# Read and parse equilibrator ΔG data. Save to Julia file. 
# """
# function mkGibbsDB(in_file_loc::String, out_file_loc::String)
#     gibbs = Dict{String, Array{Float64, 1}}()
#     GZip.open(in_file_loc) do io
#         for ln in eachline(io)
#             startswith(ln, "!") && continue
#             prts = split(ln, "\t")
#             if prts[2] != "nan"
#                 gibbs[prts[1]] = [parse(Float64, prts[2]), parse(Float64, prts[3])] 
#             end
#         end
#     end
#     JLD.save(out_file_loc, "gibbs", gibbs)
# end

# function invertdict(d)
#     nd = Dict{String, Array{String, 1}}()
#     for (k, v) in d
#         for v in d[k]
#             if haskey(nd, v)
#                 push!(nd[v], k)
#             else
#                 nd[v] = [k]
#             end
#         end
#     end
#     return nd
# end

# function add_to_dict(d, ln)
#     prts = split(ln, "\t")
#     id = string(split(prts[1], ":")[2])
#     ref = string(prts[2])
#     ref == "EMPTY" && return nothing
#     if haskey(d, id)
#         push!(d[id], ref)
#     else
#         d[id] = [ref]
#     end
# end

# function mkNameSpaceMappings(in_file_path, out_file_path)
#     bigg_to_meta = Dict{String, Array{String, 1}}()
#     kegg_to_meta = Dict{String, Array{String, 1}}()
#     open(in_file_path) do io
#         for ln in eachline(io)
#             startswith(ln, "#") && continue
#             if startswith(ln, "biggR")
#                 add_to_dict(bigg_to_meta, ln)
#             elseif startswith(ln, "keggR")
#                 add_to_dict(kegg_to_meta, ln)
#             end
#         end
#     end

#     meta_to_kegg = invertdict(kegg_to_meta)
#     bigg_to_kegg = Dict{String, Array{String, 1}}()
#     for (k, vs) in bigg_to_meta
#         for v in vs
#             if haskey(meta_to_kegg, v)
#                 if haskey(bigg_to_kegg, k)
#                     append!(bigg_to_kegg[k], meta_to_kegg[v])
#                 else
#                     bigg_to_kegg[k] = meta_to_kegg[v]
#                 end    
#             end
#         end
#     end
#     JLD.save(out_file_path, "bigg_to_kegg", bigg_to_kegg)
# end