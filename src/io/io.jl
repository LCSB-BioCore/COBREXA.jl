"""
    read_model(file_location::String))

Reads a model at `file_location` and returns a constraint based `model::CobraModel`.
Currently supported formats include SBML (.xml), Matlab (.mat) and JSON (.json) models.
The model format is inferred from the `file_location` extension.

Note, some meta-information may be lost when importing a model. Importantly, only information regarding the
reactions, metabolites and genes are imported. Currently reading JSON models captures the most meta-information
regarding reactions, metabolites and genes (e.g. the notes and annotation fields).

When importing Matlab models some annotation and notes may not be imported because of non-standard field names used by some models.
Gene reaction rules are successfully imported only if they adhere to this format: `"(YIL010W and YLR043C) or (YIL010W and YGR209C)"`,
where `or` can be interchanged with `OR, |, ||` and `and` can be interchanged with `AND, &, &&`.
Other gene reaction rules formats are not supported yet, but file an issue if your format is standard and needs to be included.

However, in all cases the basic information needed to perform constraint based analysis should be imported successfully,
e.g. stoichiometrix matrix, constraints etc..
Advanced tools that require, e.g. metabolite formulas, gene reaction rules, and KEGG or BIGG IDs, will not function if these are improperly imported.
Always inspect the imported model before running analysis (garbage in -> garbage out).
"""
function read_model(file_location::String)
    if endswith(file_location, ".json")
        try
            model = read_json_model(file_location)
        catch err
            @error "JSON model reading error.\n$err"
            model = CobraModel()
        end
    elseif endswith(file_location, ".xml")
        try
            model = reconstruct_model_sbml(file_location)
        catch err
            @error "SBML model reading error.\n$err"
            model = CobraModel()
        end
    elseif endswith(file_location, ".mat")
        try
            model = read_matlab_model(file_location)
        catch err
            @error "Matlab model reading error.\n$err"
            model = CobraModel()
        end
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
        model = CobraModel()
    end
    return model
end

"""
    save_model(model::CobraModel, file_location::String)

Save model at `file_location`. Infers format from `file_location` extension.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).

Note, only the fields contained in model are saved. Make sure that information isn't
lost between reading a model and writing a model (e.g. check gene reaction rules, notes and annotations).
"""
function save_model(model::CobraModel, file_location::String)
    if endswith(file_location, ".json")
        write_json_model(model, file_location)
    elseif endswith(file_location, ".xml")
        @warn "Not implemented!"
    elseif endswith(file_location, ".mat")
        write_matlab_model(model, file_location)
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
    end
end


"""
    parsegrr(string_rule, genes::Array{Gene, 1})

Parse a gene reaction rule string `string_rule` into a nested `gene` array `Array{Array{Gene, 1}, 1}`.

Format: (YIL010W and YLR043C) or (YIL010W and YGR209C) where `or` can also be `OR, |, ||` and where `and` can also be `AND, &, &&`.
"""
function parse_grr(s::String, genes::Array{Gene,1})
    if s == "" || isnothing(s)
        return Array{Array{Gene,1},1}()
    end
    # first get the gene id list in string format
    gene_string_rules = Array{Array{String,1},1}()
    or_genes = split(s, r"\s?(or|OR|(\|\|)|\|)\s?") # separate or terms
    for or_gene in or_genes
        and_genes = split(replace(or_gene, r"\(|\)" => ""), r"\s?(and|AND|(\&\&)|\&)\s?")
        push!(gene_string_rules, and_genes)
    end
    # now map these gene string ids to genes
    grr = Array{Array{Gene,1},1}()
    for gsr in gene_string_rules
        gene_list = Array{Gene,1}()
        for g in gsr
            gene = findfirst(genes, g)
            isnothing(gene) && (@warn "Gene not found..."; continue)
            push!(gene_list, gene)
        end
        push!(grr, gene_list)
    end
    return grr
end

"""
    unparse_grr(grr::Array{Array{Gene, 1}, 1}

Converts a nested `gene` array, `grr`, back into a grr string.
"""
function unparse_grr(grr::Array{Array{Gene,1},1})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g.id for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    return grr_string
end

"""
    reconstruct_model_sbml(file_location::String)
"""
function reconstruct_model_sbml(file_location::String)
    # m = read_sbml(file_location)
    # m is now a Model structure with:
    # m.reactions
    # m.species
    # m.compartments
    # return Model()
    return CobraModel()
end


"""
    save_sbml_model(model::CobraModel, file_location::String)
"""
function save_sbml_model(model::CobraModel, file_location::String)
    # To do...
end
