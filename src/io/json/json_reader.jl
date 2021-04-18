"""
    _read_model(file_location::String, ::Type{JSONFile})

Reads in the entire `.json` model file. 
Assumes only one model is stored in a file.
No information loss occurs, i.e. all model components can be accessed through the `JSONModel` struct, 
either generically or by explicitly accessing the internals of the struct. 

See also: [`JSONModel`](@ref)
"""
function _read_model(file_location::String, ::Type{JSONFile})
    return JSONModel(JSON.parsefile(file_location))
end
