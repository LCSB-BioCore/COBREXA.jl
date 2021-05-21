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
    workerids = [myid()],
)
    # create optimization problem, apply constraints
    save_model = :(
        begin
            optmodel = $COBREXA.make_optimization_model($model, $optimizer)
            for mod in $modifications
                mod($model, optmodel)
            end
            optmodel
        end
    )

    # load on all workers
    map(fetch, save_at.(workerids, :cobrexa_sampling_warmup_optmodel, Ref(save_model)))

    # generate warm up points, like FVA
    ret = m -> value.(m[:x]) # get all the fluxes, not just max/min like FVA
    fluxes = dpmap(
        rid -> :($COBREXA._FVA_optimize_reaction(
            cobrexa_sampling_warmup_optmodel,
            $rid,
            $ret,
        )),
        CachingPool(workerids),
        [-warmup_points warmup_points],
    )

    # get constraints from one of the workers
    lbs, ubs = get_val_from(
        workerids[1],
        :($COBREXA.get_bound_vectors(cobrexa_sampling_warmup_optmodel)),
    )

    # free the data on workers
    map(fetch, remove_from.(workerids, :cobrexa_sampling_warmup_optmodel))

    return fluxes, lbs, ubs
end
