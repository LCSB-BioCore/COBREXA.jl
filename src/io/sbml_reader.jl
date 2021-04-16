function _read_model(file_location::String, ::Type{SBMLFile}, ::Type{SBMLModel})
    return SBMLModel(SBML.readSBML(file_location)) # direct import
end

function _read_model(file_location::String, ::Type{SBMLFile}, ::Type{CoreModel})
    sbmlmodel = _read_model(file_location, SBMLFile, SBMLModel) # use other import
    return convert(CoreModel, sbmlmodel)
end
