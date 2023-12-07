
"""
    $(TYPEDSIGNATURES)

Load a FBC model representation while guessing the correct model type. Uses
`AbstractFBCModels.load`.

This overload almost always involves a search over types; do not use it in
environments where performance is critical.
"""
function load_model(path::String) where {A<:AbstractFBCModel}
    A.load(path)
end

"""
    $(TYPEDSIGNATURES)

Load a FBC model representation. Uses `AbstractFBCModels.load`.
"""
function load_model(model_type::Type{A}, path::String) where {A<:AbstractFBCModel}
    A.load(model_type, path)
end

export load_model

"""
    $(TYPEDSIGNATURES)

Save a FBC model representation. Uses `AbstractFBCModels.save`.
"""
function save_model(model::A, path::String) where {A<:AbstractFBCModel}
    A.save(model, path)
end

export save_model

"""
    $(TYPEDSIGNATURES)

Safely download a model with a known hash. All arguments are forwarded to
`AbstractFBCModels.download_data_file` -- see the documentation in the
AbstractFBCModels package for details.
"""
download_model(args...; kwargs...) = A.download_data_file(args...; kwargs...)

export download_model
