"""
    warmup_from_variability(
        n_points::Int,
        model::MetabolicModel,
        optimizer;
        modifications = [],
        workers = [myid()],
    )

Generates FVA-like warmup points for samplers, by selecting random points by
minimizing and maximizing reactions. Can not return more than 2 times the
number of reactions in the model.
"""
function warmup_from_variability(n_points::Int, model::MetabolicModel, optimizer; kwargs...)
    nr = n_reactions(model)

    n_points > 2 * nr && throw(
        DomainError(
            n_points,
            "Variability method can not generate more than $(2*nr) points from this model",
        ),
    )

    sample = shuffle(vcat(1:nr, -(1:nr)))[begin:n_points]
    warmup_from_variability(
        -filter(x -> x < 0, sample),
        filter(x -> x > 0, sample),
        model,
        optimizer;
        kwargs...,
    )
end

"""
    warmup_from_variability(
        min_reactions::Vector{Int},
        max_reactions::Vector{Int},
        model::MetabolicModel,
        optimizer;
        modifications=[],
        workers::Vector{Int}=[myid()]
    )::Matrix{Float64}

Generate FVA-like warmup points for samplers, by minimizing and maximizing the
specified reactions. The result is returned as a matrix, each point occupies as
single column in the result.
"""
function warmup_from_variability(
    min_reactions::Vector{Int},
    max_reactions::Vector{Int},
    model::MetabolicModel,
    optimizer;
    modifications = [],
    workers::Vector{Int} = [myid()],
)::Matrix{Float64}

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

    fluxes = hcat(
        dpmap(
            rid -> :($COBREXA._FVA_optimize_reaction(
                cobrexa_sampling_warmup_optmodel,
                $rid,
                optmodel -> $COBREXA.JuMP.value.(optmodel[:x]),
            )),
            CachingPool(workers),
            vcat(-min_reactions, max_reactions),
        )...,
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_sampling_warmup_optmodel))

    return fluxes
end
