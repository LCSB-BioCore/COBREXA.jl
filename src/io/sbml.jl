
"""
Import a LinearModel from a SBML file
"""
function loadSBMLModel(filename::String)
    model = SBML.readSBML(filename)

    mets, rxns, S = SBML.getS(model)
    b = zeros(length(mets))
    c = SBML.getOCs(model)

    lbu = SBML.getLBs(model)
    ubu = SBML.getUBs(model)

    # SBML has units that we can't easily interpret yet, BUT the chances are
    # that in a reasonable SBML file all units will be the same and we don't
    # need to care about them at all. So we just check it here and throw an
    # error if not.

    unit = lbu[1][2]
    getvalue = (val,_)::Tuple -> val
    getunit = (_, unit)::Tuple -> unit
    
    if any(getunit.(lbu) .!== unit) || any(getunit.(ubu) .!== unit)
        @error "The SBML file uses multiple units; loading needs conversion"
    end

    return LinearModel(S, b, c, getvalue.(lbu), getvalue.(ubu), rxns, mets)
end
