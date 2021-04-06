"""
    write_matlab_model(model::StandardModel, file_location::String)

Some information is lost here, e.g. notes and some annotations.
"""
function write_matlab_model(model::StandardModel, file_location::String)
    S = stoichiometry(model)
    b = balance(model)
    lbs, ubs = bounds(model)

    mdict = Dict(
        "c" => [r.objective_coefficient for r in model.reactions],
        "modelName" => model.id,
        "mets" => [m.id for m in model.metabolites],
        "subSystems" => [r.subsystem for r in model.reactions],
        "b" => Array(b),
        "metFormulas" => [m.formula for m in model.metabolites],
        "ub" => Array(ubs),
        "rxnNames" => [r.name for r in model.reactions],
        "description" => model.id,
        "genes" => [g.id for g in model.genes],
        "grRules" => [unparse_grr(r.grr) for r in model.reactions],
        "S" => Array(S),
        "metNames" => [m.name for m in model.metabolites],
        "lb" => Array(lbs),
        "metCharge" => [m.charge for m in model.metabolites],
        "rxns" => [r.id for r in model.reactions],
        "rxnKEGGID" => [
            join(get(r.annotation, "kegg.reaction", [""]), "; ") for r in model.reactions
        ],
        "rxnECNumbers" =>
            [join(get(r.annotation, "ec-code", [""]), "; ") for r in model.reactions],
        "rxnBiGGID" => [
            join(get(r.annotation, "bigg.reaction", [""]), "; ") for r in model.reactions
        ],
        "rxnSBOTerms" => [get(r.annotation, "sbo", "") for r in model.reactions],
        "metBiGGID" => [
            join(get(m.annotation, "bigg.metabolite", [""]), "; ") for
            m in model.metabolites
        ],
        "metSBOTerms" => [get(m.annotation, "sbo", "") for m in model.metabolites],
        "metKEGGID" => [
            join(get(m.annotation, "kegg.compound", [""]), "; ") for m in model.metabolites
        ],
        "metMetaNetXID" => [
            join(get(m.annotation, "metanetx.chemical", [""]), "; ") for
            m in model.metabolites
        ],
        "metChEBIID" =>
            [join(get(m.annotation, "chebi", [""]), "; ") for m in model.metabolites],
    )

    matwrite(file_location, Dict("model" => mdict))
end
