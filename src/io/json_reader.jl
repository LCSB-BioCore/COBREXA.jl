
"""
    load_json_model(filename::String)::JSONModel

Load and return a JSON-formatted model that is stored in `file_name`.
"""
function load_json_model(filename::String)::JSONModel
    return JSONModel(JSON.parsefile(filename))
end
