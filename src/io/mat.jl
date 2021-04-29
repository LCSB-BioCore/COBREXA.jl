
"""
    load_mat_model(file_name::String)

Load and return a MATLAB file `file_name` that contains a COBRA-compatible
model.
"""
function load_mat_model(file_name::String)::MATModel
    model_pair = first(matread(file_name))
    @_io_log @info "Loading MAT: taking a model with ID $(model_pair.first)"
    return MATModel(model_pair.second)
end

"""
    save_mat_model(model::MetabolicModel, file_name::String; model_name::String="model")

Save a [`MATModel`](@ref) in `model` to a MATLAB file `file_name` in a format
compatible with other MATLAB-based COBRA software.

In case the `model` is not `MATModel`, it will be converted automatically.

`model_name` is the identifier name for the whole model written to the MATLAB
file; defaults to just "model".
"""
function save_mat_model(model::MetabolicModel, file_path::String; model_name = "model")
    m =
        typeof(model) == MATModel ? model :
        begin
            @_io_log @warn "Automatically converting $(typeof(model)) to MATModel for saving, information may be lost."
            convert(MATModel, model)
        end
    matwrite(file_path, Dict(model_name => m.mat))
end
