"""
    warmup_from_variability(
        n_points::Int,
        model::MetabolicModel,
        optimizer;
        modifications = [],
        workers = [myid()],
    )

Generates warmup points for samplers by minimizing and maximizing random
reactions. Very similar to
[`flux_variability_analysis`](@ref).
"""
function warmup_from_variability(
    n_points::Int,
    model::MetabolicModel,
    optimizer;
    kwargs...
)

    
end

function warmup_from_variability(
    min_reactions::Vector{Int},
    max_reactions::Vector{Int},
    model::MetabolicModel,
    optimizer;
    modifications=[],
    workers::Vector{Int}=[myid()]
)::Tuple{Matrix{Float64}, AbstractVector{Float64}, AbstractVector{Float64}}
    # create optimization problem at workers, apply modifications
    save_model = :(
        begin
            model = $model
            optmodel = $COBREXA.make_optimization_model(model, $optimizer)
            for mod in $modifications
                mod(model, optmodel)
            end
            optmodel
        end
    )

    map(fetch, save_at.(workers, :cobrexa_sampling_warmup_optmodel, Ref(save_model)))

    fluxes = dpmap(
        rid -> :($COBREXA._FVA_optimize_reaction(
            cobrexa_sampling_warmup_optmodel,
            $rid,
            optmodel -> value.(optmodel[:x]),
        )),
        CachingPool(workers),
        vcat(-min_reactions, max_reactions),
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_sampling_warmup_optmodel))

    return hcat(fluxes...)
end
