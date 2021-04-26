
"""
    load_sbml_model(file_name::String)::SBMLModel

Load and return a SBML XML model in `file_name`.
"""
function load_sbml_model(file_name::String)::SBMLModel
    return SBMLModel(SBML.readSBML(file_name))
end
