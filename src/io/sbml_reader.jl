function _read_sbml_model(file_location::String, ::SBMLModel)
    SBMLModel(SBML.readSBML(file_location))
end

function _read_sbml_model(file_location::String, ::StandardModel)
    # m = read_sbml(file_location)
    # m is now a Model structure with:
    # m.reactions
    # m.species
    # m.compartments
    # return Model()
    return StandardModel()
end
