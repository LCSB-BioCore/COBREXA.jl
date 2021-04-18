"""
    _read_model(file_location::String, ::Type{MATFile})

Reads in the entire `.mat` model file. 
Assumes only one model is stored in a file.
No information loss occurs, i.e. all model components can be accessed through the `MATModel` struct, 
either generically or by explicitly accessing the internals of the struct. 

See also: [`MATModel`](@ref)
"""
function _read_model(file_location::String, ::Type{MATFile})
    matfile = matread(file_location)
    model_name = collect(keys(matfile))[1] # assume only one model per m-file
    return MATModel(model_name, matfile[model_name])
end
