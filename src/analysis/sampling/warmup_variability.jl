"""
$(TYPEDSIGNATURES)

Generates FVA-like warmup points for samplers, by selecting random points by
minimizing and maximizing reactions. Can not return more than 2 times the
number of reactions in the model.
"""
function warmup_from_variability(
    model::AbstractMetabolicModel,
    optimizer,
    n_points::Int,
    seed = rand(Int);
    kwargs...,
)
    nr = n_reactions(model)

    n_points > 2 * nr && throw(
        DomainError(
            n_points,
            "Variability method can not generate more than $(2*nr) points from this model",
        ),
    )

    sample = shuffle(StableRNG(seed % UInt), vcat(1:nr, -(1:nr)))[begin:n_points]
    warmup_from_variability(
        model,
        optimizer,
        -filter(x -> x < 0, sample),
        filter(x -> x > 0, sample);
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

Generate FVA-like warmup points for samplers, by minimizing and maximizing the
specified reactions. The result is returned as a matrix, each point occupies as
single column in the result.

!!! warning Limited effect of modifications in `warmup_from_variability`
    Modifications of the optimization model applied in `modifications`
    parameter that change the semantics of the model have an effect on the
    warmup points, but do not automatically carry to the subsequent sampling.
    Users are expected to manually transplant any semantic changes to the
    actual sampling functions, such as [`affine_hit_and_run`](@ref).
"""
function warmup_from_variability(
    model::AbstractMetabolicModel,
    optimizer,
    min_reactions::AbstractVector{Int} = 1:n_reactions(model),
    max_reactions::AbstractVector{Int} = 1:n_reactions(model);
    modifications = [],
    workers::Vector{Int} = [myid()],
)::Matrix{Float64}

    # create optimization problem at workers, apply modifications
    save_model = :(
        begin
            local model = $model
            local optmodel = $make_optimization_model(model, $optimizer)
            for mod in $modifications
                mod(model, optmodel)
            end
            optmodel
        end
    )

    asyncmap(fetch, save_at.(workers, :cobrexa_sampling_warmup_optmodel, Ref(save_model)))

    fluxes = hcat(
        dpmap(
            rid -> :($_maximize_warmup_reaction(
                cobrexa_sampling_warmup_optmodel,
                $rid,
                om -> $JuMP.value.(om[:x]),
            )),
            CachingPool(workers),
            vcat(-min_reactions, max_reactions),
        )...,
    )

    # free the data on workers
    asyncmap(fetch, remove_from.(workers, :cobrexa_sampling_warmup_optmodel))

    return fluxes
end

"""
$(TYPEDSIGNATURES)

A helper function for finding warmup points from reaction variability.
"""
function _maximize_warmup_reaction(opt_model, rid, ret)
    sense = rid > 0 ? MAX_SENSE : MIN_SENSE
    var = all_variables(opt_model)[abs(rid)]

    @objective(opt_model, sense, var)
    optimize!(opt_model)

    is_solved(opt_model) ? ret(opt_model) : nothing
end
