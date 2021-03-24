
"""
Import a LinearModel from a SBML file
"""
loadSBMLModel(filename::String)::SBMLModel = SBMLModel(SBML.readSBML(filename))
