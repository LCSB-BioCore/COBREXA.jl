
"""
$(TYPEDSIGNATURES)

Load and return a MATLAB file `file_name` that contains a COBRA-compatible
model.
"""
function load_mat_model(file_name::String)::MATModel
    model_pair = first(matread(file_name))
    @io_log @info "Loading MAT: taking a model with ID $(model_pair.first)"
    return MATModel(model_pair.second)
end

"""
$(TYPEDSIGNATURES)

Save a [`MATModel`](@ref) in `model` to a MATLAB file `file_name` in a format
compatible with other MATLAB-based COBRA software.

In case the `model` is not `MATModel`, it will be converted automatically.

`model_name` is the identifier name for the whole model written to the MATLAB
file; defaults to just "model".
"""
function save_mat_model(model::AbstractMetabolicModel, file_path::String; model_name = "model")
    m =
        typeof(model) == MATModel ? model :
        begin
            @io_log @warn "Automatically converting $(typeof(model)) to MATModel for saving, information may be lost."
            convert(MATModel, model)
        end
    matwrite(file_path, Dict(model_name => m.mat))
end
