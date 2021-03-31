
"""
Import a LinearModel from a SBML file
"""
load_sbml_model(filename::String)::SBMLModel = SBMLModel(SBML.readSBML(filename))
