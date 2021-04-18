"""
    _read_model(file_location::String, ::Type{SBMLFile})

Reads in the entire `.xml` model file. 
Assumes only one model is stored in a file.
No information loss occurs, i.e. all model components can be accessed through the `SBMLModel` struct, 
either generically or by explicitly accessing the internals of the struct. 

See also: [`SBMLModel`](@ref)
"""
function _read_model(file_location::String, ::Type{SBMLFile})
    return SBMLModel(SBML.readSBML(file_location)) # direct import
end
