"""
$(TYPEDSIGNATURES)

Generic function for loading models that chooses a specific loader function
based on the `extension` argument (e.g., `".xml"` chooses loading of the SBML
model format), or throws an error.  By default the extension from `file_name`
is used.

Currently, these model types are supported:

- SBML models (`*.xml`, loaded with [`load_sbml_model`](@ref))
- JSON models (`*.json`, loaded with [`load_json_model`](@ref))
- MATLAB models (`*.mat`, loaded with [`load_mat_model`](@ref))
- HDF5 models (`*.h5`, loaded with [`load_h5_model`](@ref))
"""
function load_model(
    file_name::String;
    extension = last(splitext(file_name)),
)::MetabolicModel

    if extension == ".json"
        return load_json_model(file_name)
    elseif extension == ".xml"
        return load_sbml_model(file_name)
    elseif extension == ".mat"
        return load_mat_model(file_name)
    elseif extension == ".h5"
        return load_h5_model(file_name)
    else
        throw(DomainError(extension, "Unknown file extension"))
    end
end

"""
$(TYPEDSIGNATURES)

Helper function that loads the model using [`load_model`](@ref) and returns it
converted to `type`.

# Example:

    load_model(CoreModel, "mySBMLModel.xml")
"""
function load_model(
    type::Type{T},
    file_name::String;
    extension = last(splitext(file_name)),
)::T where {T<:MetabolicModel}
    convert(type, load_model(file_name; extension))
end

"""
$(TYPEDSIGNATURES)

Generic function for saving models that chooses a specific writer function from
the `extension` argument (such as `".xml"` for SBML format), or throws an
error. By default the extension from `file_name` is used.

Currently, these model types are supported:

- SBML models (`*.xml`, saved with [`save_sbml_model`](@ref))
- JSON models (`*.json`, saved with [`save_json_model`](@ref))
- MATLAB models (`*.mat`, saved with [`save_mat_model`](@ref))
- HDF5 models (`*.h5`, saved with [`save_h5_model`](@ref))
"""
function save_model(
    model::MetabolicModel,
    file_name::String;
    extension = last(splitext(file_name)),
)
    if extension == ".json"
        return save_json_model(model, file_name)
    elseif extension == ".xml"
        return save_sbml_model(model, file_name)
    elseif extension == ".mat"
        return save_mat_model(model, file_name)
    elseif extension == ".h5"
        return save_h5_model(model, file_name)
    else
        throw(DomainError(extension, "Unknown file extension"))
    end
end
