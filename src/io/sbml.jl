
"""
    load_sbml_model(file_name::String)::SBMLModel

Load and return a SBML XML model in `file_name`.
"""
function load_sbml_model(file_name::String)::SBMLModel
    return SBMLModel(SBML.readSBML(file_name))
end

"""
    write_sbml_model(model::MetabolicModel, file_name::String)

Write a given SBML model to `file_name`.
"""
function save_sbml_model(model::MetabolicModel, file_name::String)
    m =
        typeof(model) == SBMLModel ? model :
        begin
            @_io_log @warn "Automatically converting $(typeof(model)) to SBMLModel for saving, information may be lost."
            convert(SBMLModel, model)
        end

    SBML.writeSBML(m.sbml, file_name)
end
