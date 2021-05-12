"""
warmup(
    model::MetabolicModel,
    optimizer;
    warmup_points::Vector{Int},
    modifications = [],
    workers = [myid()],

Generates warmup points for samplers by sequentially minimizing and maximizing 
reactions at `warmup_points`. Very similar to [`flux_variability_analysis`](@ref).
"""
function warmup(
    model::MetabolicModel,
    optimizer;
    warmup_points::Vector{Int},
    modifications = [],
    workers = [myid()],
)
    # create optimization problem, apply constraints, load on all workers
    save_model = :(
        begin
            optmodel = $COBREXA.make_optimization_model($model, $optimizer)
            for mod in $modifications
                mod($model, optmodel)
            end
            optmodel
        end
    )

    map(fetch, save_at.(workers, :cobrexa_hit_and_run_warmup_model, Ref(save_model)))

    ret = m -> value.(m[:x]) # get all the fluxes
    fluxes = dpmap(
        rid -> :($COBREXA._FVA_optimize_reaction(
            cobrexa_hit_and_run_warmup_model,
            $rid,
            $ret,
        )),
        CachingPool(workers),
        [-warmup_points warmup_points],
    )

    # get constraints from one of the workers
    lbs, ubs = get_val_from(
        workers[1],
        :(COBREXA.get_bound_vectors(cobrexa_hit_and_run_warmup_model)),
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_hit_and_run_warmup_model))

    return fluxes, lbs, ubs
end
