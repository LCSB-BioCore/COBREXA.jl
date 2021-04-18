
"""
    _infer_best_model_type(filename::String)::Maybe{Type}

Guess the type of the model that's best for representing stuff in a file named
`file_name`.
"""
function _infer_best_model_type(filename::String)::Maybe{Type}
    if endswith(file_name, ".json")
        return StandardModel
    elseif endswith(file_name, ".xml")
        return SBMLModel
    elseif endswith(file_name, ".mat")
        return MATModel
    else
        return nothing
    end
end

"""
    read_model(file_location::String, ::Type{StandardModel})

Reads a model at `file_location` and returns a constraint based `model::StandardModel`.
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
function read_model(filename::String, type::Maybe{Type} = nothing)
    if isnothing(type)
        type = _infer_best_model_type(filename)
    end
    if isnothing(type) # still!
        throw(DomainError(filename, "Could not infer a proper model type to load from file extension"))
    end
    return read_model(filename, type)
end

"""
    write_model(model::StandardModel, file_location::String)

Save model at `file_location`. Infers format from `file_location` extension.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).

Note, only the fields contained in model are saved. Make sure that information isn't
lost between reading a model and writing a model (e.g. check gene reaction rules, notes and annotations).
"""
function write_model(model::MetabolicModel, file_location::String)
    if isnothing(type)
        type = _infer_best_model_type
    end
    if isnothing(type) # still!
        throw(DomainError(filename, "Could not infer a proper model type to write from file extension"))
    end
    if inferred_type == UNKNOWNFile
        @warn "File type not supported."
        return nothing
    else
        _write_model(model, inferred_type, file_location)
        return nothing
    end
end
