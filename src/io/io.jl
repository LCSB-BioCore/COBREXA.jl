
"""
$(TYPEDSIGNATURES)

Generic function for loading models that chooses a specific loader function
from the `file_name` extension, or throws an error.

Currently, these model types are supported:

- SBML models (`*.xml`, loaded with [`load_sbml_model`](@ref))
- JSON models (`*.json`, loaded with [`load_json_model`](@ref))
- MATLAB models (`*.mat`, loaded with [`load_mat_model`](@ref))
- HDF5 models (`*.h5`, loaded with [`load_h5_model`](@ref))
"""
function load_model(file_name::String)::AbstractMetabolicModel
    if endswith(file_name, ".json")
        return load_json_model(file_name)
    elseif endswith(file_name, ".xml")
        return load_sbml_model(file_name)
    elseif endswith(file_name, ".mat")
        return load_mat_model(file_name)
    elseif endswith(file_name, ".h5")
        return load_h5_model(file_name)
    else
        throw(DomainError(file_name, "Unknown file extension"))
    end
end


"""
$(TYPEDSIGNATURES)

Helper function tht loads the model using [`load_model`](@ref) and return it
converted to `type`.

# Example:

    load_model(CoreModel, "mySBMLModel.xml")
"""
function load_model(type::Type{T}, file_name::String)::T where {T<:AbstractMetabolicModel}
    convert(type, load_model(file_name))
end

"""
$(TYPEDSIGNATURES)

Generic function for saving models that chooses a specific writer function
from the `file_name` extension, or throws an error.

Currently, these model types are supported:

- SBML models (`*.xml`, saved with [`save_sbml_model`](@ref))
- JSON models (`*.json`, saved with [`save_json_model`](@ref))
- MATLAB models (`*.mat`, saved with [`save_mat_model`](@ref))
- HDF5 models (`*.h5`, saved with [`save_h5_model`](@ref))
"""
function save_model(model::AbstractMetabolicModel, file_name::String)
    if endswith(file_name, ".json")
        return save_json_model(model, file_name)
    elseif endswith(file_name, ".xml")
        return save_sbml_model(model, file_name)
    elseif endswith(file_name, ".mat")
        return save_mat_model(model, file_name)
    elseif endswith(file_name, ".h5")
        return save_h5_model(model, file_name)
    else
        throw(DomainError(file_name, "Unknown file extension"))
    end
end
