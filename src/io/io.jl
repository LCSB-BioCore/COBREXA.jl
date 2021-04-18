struct JSONFile end
struct SBMLFile end
struct MATFile end
struct UNKNOWNFile end

"""
    _infer_file_type(file_name::String)

Infer the file type given its extension.
"""
function _infer_file_type(file_name::String)
    if endswith(file_name, ".json")
        return JSONFile
    elseif endswith(file_name, ".xml")
        return SBMLFile
    elseif endswith(file_name, ".mat") 
        return MATFile
    end
    return UNKNOWNFile
end

"""
    read_model(file_location::String)

Read a model from file at `file_location`. 
This function infers the model type based on the file extension.
Files ending with `.xml` are read as `SBMLModel`s, files ending with `.mat` are read as `MATModel`s, 
and files ending with `.json` are read as `JSONModel`s.

Note, that all analysis functions work with these model types, although high performance computations
should be performed by converting these models to either `StandardModel`, `CoreModel` or `CoreModelCoupled`. 
Also note, no information loss occurs when importing models, but information loss may occur when converting
between different model types (e.g. non-conventional meta-information).
"""
function read_model(file_location::String)
    inferred_type = _infer_file_type(file_location)
    if inferred_type == UNKNOWNFile
        @warn "File type not supported."
        return nothing
    else
        return _read_model(file_location, inferred_type)
    end
end

"""
    write_model(model::StandardModel, file_location::String)

Save model at `file_location`. 
Infers format from `file_location` extension.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).

Note, only the fields contained in model are saved. Make sure that information isn't
lost between reading a model and writing a model (e.g. check gene reaction rules, notes and annotations).
"""
function write_model(model::MetabolicModel, file_location::String)
    inferred_type = _infer_file_type(file_location)
    if inferred_type == UNKNOWNFile
        @warn "File type not supported."
        return nothing
    else
        _write_model(model, inferred_type, file_location)
        return nothing
    end
end
