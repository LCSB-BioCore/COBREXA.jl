
@testset "Conversion from and to SBML model" begin
    sbmlm = load_sbml_model(model_paths["ecoli_core_model.xml"])
    sm = convert(StandardModel, sbmlm)
    sbmlm2 = convert(SBMLModel, sm)

    @test Set(reactions(sbmlm)) == Set(reactions(sbmlm2))
    @test Set(reactions(sbmlm)) == Set(reactions(sm))
    @test Set(metabolites(sbmlm)) == Set(metabolites(sbmlm2))
    @test all([
        sbmlm.sbml.reactions[i].stoichiometry == sbmlm2.sbml.reactions[i].stoichiometry for
        i in reactions(sbmlm2)
    ])
end
