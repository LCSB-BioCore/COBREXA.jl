
@testset "Conversion from and to SBML model" begin
    sbmlm = load_sbml_model(model_paths["ecoli_core_model.xml"])
    sm = convert(ObjectModel, sbmlm)
    sbmlm2 = convert(SBMLModel, sm)

    @test Set(variables(sbmlm)) == Set(variables(sbmlm2))
    @test Set(variables(sbmlm)) == Set(variables(sm))
    @test Set(metabolites(sbmlm)) == Set(metabolites(sbmlm2))
    sp(x) = x.species
    @test all([
        issetequal(
            sp.(sbmlm.sbml.reactions[i].reactants),
            sp.(sbmlm2.sbml.reactions[i].reactants),
        ) && issetequal(
            sp.(sbmlm.sbml.reactions[i].products),
            sp.(sbmlm2.sbml.reactions[i].products),
        ) for i in variables(sbmlm2)
    ])
    st(x) = isnothing(x.stoichiometry) ? 1.0 : x.stoichiometry
    @test all([
        issetequal(
            st.(sbmlm.sbml.reactions[i].reactants),
            st.(sbmlm2.sbml.reactions[i].reactants),
        ) && issetequal(
            st.(sbmlm.sbml.reactions[i].products),
            st.(sbmlm2.sbml.reactions[i].products),
        ) for i in variables(sbmlm2)
    ])
end

@testset "SBMLModel generic interface" begin
    model = load_model(model_paths["e_coli_core.xml"])

    @test reaction_stoichiometry(model, "R_EX_ac_e") == Dict("M_ac_e" => -1)
end