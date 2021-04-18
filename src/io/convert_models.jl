"""
Convert a model into a `CoreModel`
"""
function Base.convert(::Type{CoreModel}, model::MetabolicModel)
    lb, ub = bounds(model)
    return CoreModel(stoichiometry(model), balance(model), c=objective(model), lb, ub, reactions(model), metabolites(model))
end

function Base.convert(::Type{StandardModel}, model::JSONModel)
# function _read_model(file_location::String, ::Type{JSONFile}, ::Type{StandardModel})
#     modeldict = JSON.parsefile(file_location)

#     modelid = modeldict["id"]

#     mets = Metabolite[]
#     for i in modeldict["metabolites"]
#         met = Metabolite()

#         for (k, v) in i
#             if k == "id"
#                 met.id = v
#             elseif k == "name"
#                 met.name = v
#             elseif k == "formula"
#                 met.formula = v
#             elseif k == "charge"
#                 met.charge = v
#             elseif k == "compartment"
#                 met.compartment = v
#             elseif k == "notes"
#                 met.notes = Dict{String,Vector{String}}(kk => vv for (kk, vv) in v)
#             elseif k == "annotation"
#                 for (kk, vv) in v
#                     if typeof(vv) == String
#                         met.annotation[kk] = vv
#                     else
#                         met.annotation[kk] = convert(Vector{String}, vv)
#                     end
#                 end
#             else
#                 @warn "Unrecognized metabolite field: $k"
#             end
#         end
#         push!(mets, met)
#     end

#     genes = Gene[]
#     for i in modeldict["genes"]
#         gene = Gene()
#         for (k, v) in i
#             if k == "id"
#                 gene.id = v
#             elseif k == "name"
#                 gene.name = v
#             elseif k == "notes"
#                 for (kk, vv) in v
#                     gene.notes[kk] = vv
#                 end
#             elseif k == "annotation"
#                 for (kk, vv) in v
#                     if typeof(vv) == String
#                         gene.annotation[kk] = vv
#                     else
#                         gene.annotation[kk] = convert(Vector{String}, vv)
#                     end
#                 end
#             else
#                 @warn "Unrecognized gene field: $k"
#             end
#         end
#         push!(genes, gene)
#     end

#     rxns = Reaction[]
#     for i in modeldict["reactions"]
#         rxn = Reaction()
#         for (k, v) in i
#             if k == "id"
#                 rxn.id = v
#             elseif k == "name"
#                 rxn.name = v
#             elseif k == "metabolites"
#                 for (kk, vv) in v
#                     ind = findfirst(x -> x.id == kk, mets)
#                     if isnothing(ind)
#                         @warn "Metabolite $kk not found in reaction assignment."
#                         continue
#                     else
#                         rxn.metabolites[mets[ind]] = vv
#                     end
#                 end
#             elseif k == "lower_bound"
#                 rxn.lb = v
#             elseif k == "upper_bound"
#                 rxn.ub = v
#             elseif k == "gene_reaction_rule"
#                 rxn.grr = _parse_grr(v, genes)
#             elseif k == "subsystem"
#                 rxn.subsystem = v
#             elseif k == "notes"
#                 for (kk, vv) in v
#                     rxn.notes[kk] = vv
#                 end
#             elseif k == "annotation"
#                 for (kk, vv) in v
#                     if typeof(vv) == String
#                         rxn.annotation[kk] = vv
#                     else
#                         rxn.annotation[kk] = convert(Vector{String}, vv)
#                     end
#                 end
#             elseif k == "objective_coefficient"
#                 rxn.objective_coefficient = v
#             else
#                 @warn "Unrecognized reaction field: $k"
#             end
#         end
#         push!(rxns, rxn)
#     end

#     return StandardModel(modelid, rxns, mets, genes)
# end
end

function Base.convert(::Type{StandardModel}, model::MATModel)
# function _read_model(file_location::String, ::Type{MFile}, ::Type{StandardModel})
#     matfile = matread(file_location)
#     model_name = collect(keys(matfile))[1]
#     modeldict = matfile[model_name]

#     # the model_id can be written in many places, try varying levels of specificity
#     model_id = haskey(modeldict, "description") ? modeldict["description"] : model_name
#     model_id = haskey(modeldict, "modelName") ? modeldict["modelName"] : model_name # more specific

#     mets = Metabolite[]
#     for i in eachindex(modeldict["mets"])
#         met = Metabolite()

#         if haskey(modeldict, "mets")
#             met.id = modeldict["mets"][i]
#         else
#             continue
#         end

#         if haskey(modeldict, "metNames")
#             met.name = modeldict["metNames"][i]
#         end

#         if haskey(modeldict, "metFormulas")
#             met.formula = string(modeldict["metFormulas"][i])
#         end

#         if haskey(modeldict, "metFormula")
#             met.formula = string(modeldict["metFormula"][i])
#         end

#         if haskey(modeldict, "metCharge") && !isnan(modeldict["metCharge"][i])
#             met.charge = modeldict["metCharge"][i]
#         end

#         if haskey(modeldict, "metCharges") && !isnan(modeldict["metCharges"][i])
#             met.charge = modeldict["metCharges"][i]
#         end

#         # look for annotation data, assume delimited by "; "
#         anno_kid = Dict(
#             "metBiGGID" => "bigg.metabolite",
#             "metKEGGID" => "kegg.compound",
#             "metMetaNetXID" => "metanetx.chemical",
#             "metChEBIID" => "chebi",
#         )
#         for (anno, kid) in anno_kid
#             if haskey(modeldict, anno)
#                 met.annotation[kid] = string.(split(string(modeldict[anno][i]), "; "))
#             end
#         end

#         if haskey(modeldict, "metSBOTerms")
#             met.annotation["sbo"] = string(modeldict["metSBOTerms"][i])
#         end

#         if haskey(modeldict, "metNotes")
#             met.notes["note"] = string.(split(string(modeldict["metNotes"][i]), "; "))
#         end

#         push!(mets, met)
#     end

#     genes = Gene[]
#     for i in eachindex(modeldict["genes"])
#         gene = Gene()
#         if haskey(modeldict, "genes")
#             gene.id = modeldict["genes"][i]
#         else
#             continue
#         end

#         # name and other fields don't generally exist in the matlab models,
#         # ignoring them for now

#         push!(genes, gene)
#     end

#     rxns = Reaction[]
#     for i in eachindex(modeldict["rxns"])
#         rxn = Reaction()
#         if haskey(modeldict, "rxns")
#             rxn.id = modeldict["rxns"][i]
#         else
#             continue
#         end

#         if haskey(modeldict, "rxnNames")
#             rxn.name = modeldict["rxnNames"][i]
#         end

#         metinds = findall(x -> x .!= 0.0, modeldict["S"][:, i])
#         rxn.metabolites =
#             Dict{Metabolite,Float64}(mets[j] => modeldict["S"][j, i] for j in metinds)

#         if haskey(modeldict, "lb")
#             rxn.lb = modeldict["lb"][i]
#         end

#         if haskey(modeldict, "ub")
#             rxn.ub = modeldict["ub"][i]
#         end

#         if haskey(modeldict, "grRules")
#             grr_string = modeldict["grRules"][i]
#             rxn.grr = _parse_grr(grr_string, genes)
#         end

#         rxn.subsystem = join(modeldict["subSystems"][i], "; ")

#         if haskey(modeldict, "c")
#             rxn.objective_coefficient = modeldict["c"][i]
#         end

#         # look for some annotations
#         anno_kids = Dict(
#             "rxnKEGGID" => "kegg.reaction",
#             "rxnECNumbers" => "ec-code",
#             "rxnBiGGID" => "bigg.reaction",
#         )

#         for (anno, kid) in anno_kids
#             if haskey(modeldict, anno)
#                 rxn.annotation[kid] = string.(split(string(modeldict[anno][i]), "; "))
#             end
#         end

#         if haskey(modeldict, "rxnSBOTerms")
#             rxn.annotation["sbo"] = string(modeldict["rxnSBOTerms"][i])
#         end

#         # look for some notes
#         if haskey(modeldict, "rxnNotes")
#             rxn.notes["note"] = string.(split(string(modeldict["rxnNotes"][i]), "; "))
#         end

#         push!(rxns, rxn)
#     end

#     return StandardModel(model_id, rxns, mets, genes)
# end

end

function Base.convert(::Type{StandardModel}, model::SBMLModel)

end

# """
# Convert a CoreModel to exportable format
# SparseVectors are not written and read properly, SparseMatrix is okay
# """
# function Base.convert(::AbstractDict, model::CoreModel)
#     xl, xu = bounds(model)
#     return Dict(
#         "S" => stoichiometry(model),
#         "b" => Vector(balance(model)),
#         "c" => Vector(objective(model)),
#         "ub" => Vector(xu),
#         "lb" => Vector(xl),
#         "rxns" => reactions(model),
#         "mets" => metabolites(model),
#     )
# end
