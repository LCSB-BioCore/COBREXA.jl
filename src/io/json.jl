"""
    load_json_model(filename::String)::JSONModel

Load and return a JSON-formatted model that is stored in `file_name`.
"""
function load_json_model(filename::String)::JSONModel
    return JSONModel(JSON.parsefile(filename))
end

"""
    save_json_model(model::MetabolicModel, file_name::String)

Save a [`JSONModel`](@ref) in `model` to a JSON file `file_name`.

In case the `model` is not `JSONModel`, it will be converted automatically.
"""
function save_json_model(model::MetabolicModel, file_name::String)
    m =
        typeof(model) == JSONModel ? model :
        begin
            @_io_log @warn "Automatically converting $(typeof(model)) to JSONModel for saving, information may be lost."
            convert(JSONModel, model)
        end

    JSON.print(open(file_name, "w"), m.m)
end
