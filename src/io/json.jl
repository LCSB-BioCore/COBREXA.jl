"""
$(TYPEDSIGNATURES)

Load and return a JSON-formatted model that is stored in `file_name`.
"""
function load_json_model(filename::String)::JSONModel
    return JSONModel(JSON.parsefile(filename))
end

"""
$(TYPEDSIGNATURES)

Save a [`JSONModel`](@ref) in `model` to a JSON file `file_name`.

In case the `model` is not `JSONModel`, it will be converted automatically.
"""
function save_json_model(model::AbstractMetabolicModel, file_name::String)
    m =
        typeof(model) == JSONModel ? model :
        begin
            @io_log @warn "Automatically converting $(typeof(model)) to JSONModel for saving, information may be lost."
            convert(JSONModel, model)
        end

    open(f -> JSON.print(f, m.json), file_name, "w")
end
