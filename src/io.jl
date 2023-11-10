
"""
    $(TYPEDSIGNATURES)

Load a FBC model representation while guessing the correct model type. Uses
`AbstractFBCModels.load`.

This overload almost always involves a search over types; do not use it in
environments where performance is critical.
"""
function load_fbc_model(path::String) where A<:AbstractFBCModel
    A.load(path)
end

"""
    $(TYPEDSIGNATURES)

Load a FBC model representation. Uses `AbstractFBCModels.load`.
"""
function load_fbc_model(model_type::Type{A}, path::String) where A<:AbstractFBCModel
    A.load(model_type, path)
end

export load_fbc_model

"""
    $(TYPEDSIGNATURES)

Save a FBC model representation. Uses `AbstractFBCModels.save`.
"""
function save_fbc_model(model::A, path::String) where A<:AbstractFBCModel
    A.save(model, path)
end

export save_fbc_model
