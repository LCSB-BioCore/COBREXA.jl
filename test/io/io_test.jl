@testset "IO models from file" begin

    # E. coli models - realistic size models
    iJO1366_xml = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.xml",
        joinpath("data", "iJO1366.xml"),
        "d6d9ec61ef6f155db5bb2f49549119dc13b96f6098b403ef82ea4240b27232eb",
    )
    sbmlmodel_ecoli = read_model(iJO1366_xml)
    @test n_reactions(sbmlmodel_ecoli) == 2583

    iJO1366_mat = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        joinpath("data", "iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    )
    matlabmodel_ecoli = read_model(iJO1366_mat)
    @test n_reactions(matlabmodel_ecoli) == 2583

    iJO1366_json = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        joinpath("data", "iJO1366.json"),
        "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
    )
    jsonmodel_ecoli = read_model(iJO1366_json)
    @test n_reactions(matlabmodel_ecoli) == 2583

    # test if same reactions and metabolites are read
    @test length(intersect(reactions(matlabmodel_ecoli), reactions(jsonmodel_ecoli))) == 0
    @test length(intersect(reactions(sbmlmodel_ecoli), reactions(jsonmodel_ecoli))) == 0
    @test length(intersect(metabolites(matlabmodel_ecoli), metabolites(jsonmodel_ecoli))) == 0
    @test length(intersect(metabolites(sbmlmodel_ecoli), metabolites(jsonmodel_ecoli))) == 0
    
    # test if stoichiometric matrices are the same
    sbml_S = sum(stoichiometric(sbmlmodel_ecoli), dims=(1,2))[1]
    json_S = sum(stoichiometric(sbmlmodel_ecoli), dims=(1,2))[1]
    mat_S = sum(stoichiometric(sbmlmodel_ecoli), dims=(1,2))[1]
    @test sbml_S == json_S
    @test sbml_S == mat_S    

    # test if bounds are the same
    sbml_bounds = Dict(zip(reactions(sbmlmodel_ecoli), bounds(sbmlmodel_ecoli)))
    json_bounds = Dict(zip(reactions(sbmlmodel_ecoli), bounds(sbmlmodel_ecoli)))
    mat_bounds = Dict(zip(reactions(sbmlmodel_ecoli), bounds(sbmlmodel_ecoli)))
    bound_check = true
    for k in reactions(sbmlmodel_ecoli)
        if !all(sbml_bounds[k] .== json_bounds[k]) || !all(sbml_bounds[k] .== mat_bounds[k])
            bound_check = false
        end
    end
    @test bound_check
end
