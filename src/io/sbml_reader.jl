function _read_model(file_location::String, ::Type{SBMLFile}, ::Type{SBMLModel})
    return SBMLModel(SBML.readSBML(file_location)) # direct import
end

function _read_model(file_location::String, ::Type{SBMLFile}, ::Type{LinearModel})
    sbmlmodel = _read_model(file_location, SBMLFile, SBMLModel) # use other import
    return convert(LinearModel, sbmlmodel)
end

function _read_model(file_location::String, ::Type{SBMLFile}, ::Type{StandardModel})
    m = _read_model(file_location, SBMLFile, SBMLModel) # use other import
    # do more processing
    
    return StandardModel()
end
