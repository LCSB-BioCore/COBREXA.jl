
"""
    screen_model_variants(
        model::MetabolicModel,
        model_variants::Array,
        analysis;
        workers = [myid()],
    )::Array

Take vectors of model modifications in `model_variants` array and execute the
function `analysis` on all variants of the `model` specified by
`model_variants`. The computation is distributed over worker IDs in `workers`.

The array of variants must contain vectors of single-parameter functions, these
are applied to model in order. The functions must *not* modify the model, but
rather return a modified copy. The copy should be made as shallow as possible,
to increase memory efficiency of the process. Modifications that modify the
argument model in-place will cause unpredictable results. Refer to the
definition of [`screen_one_variant`](@ref) for details.

The function `analysis` will receive a single argument (the modified model).
The results of running `analysis` are collected in to the resulting array, in a
way that preserves the shape of the `model_variants`, similarly as with `pmap`.

The modification and analysis functions are transferred to `workers` as-is; all
packages required to run them (e.g. the optimization solvers) must be loaded
there. Typically, you want to use the macro `@everywhere using
MyFavoriteSolver` from `Distributed` package for loading the solvers.

# Example
```
function reverse_reaction(i::Int)
    (model::CoreModel) -> begin
        mod = copy(model)
        mod.S[:,i] .*= -1   # this is unrealistic but sufficient for demonstration
        mod
    end
end

m = load_model(CoreModel, "e_coli_core.xml")

screen_model_variants(m,
           [
               [reverse_reaction(5)],
               [reverse_reaction(3), reverse_reaction(6)]
           ],
           mod -> mod.S[:,3])  # observe the changes in S

screen_model_variants(m,
    [
        [reverse_reaction(5)],
        [reverse_reaction(3), reverse_reaction(6)]
    ],
    mod -> flux_balance_analysis_vec(mod, GLPK.Optimizer))  # run analysis
"""
function screen_model_variants(
    model::MetabolicModel,
    model_variants::Array,
    analysis;
    workers = [myid()],
)::Array

    map(fetch, save_at.(workers, :cobrexa_screen_model, Ref(model)))
    map(fetch, save_at.(workers, :cobrexa_screen_analysis_fn, Ref(analysis)))
    map(fetch, get_from.(workers, Ref(:(precache!(cobrexa_screen_model)))))

    res = dpmap(
        mods -> :($COBREXA.screen_one_variant(
            cobrexa_screen_model,
            $mods,
            cobrexa_screen_analysis_fn,
        )),
        CachingPool(workers),
        model_variants,
    )

    map(fetch, remove_from.(workers, :cobrexa_screen_model))
    map(fetch, remove_from.(workers, :cobrexa_screen_analysis_fn))

    return res
end

"""
    screen_one_variant(model::MetabolicModel, model_modifications, analysis)

Helper function for [`screen_model_variants`](@ref) applies all single-argument
functions in `model_modifications` to the `model` (in order from "first" to
"last"), and executes `analysis` on the result.
"""
function screen_one_variant(model::MetabolicModel, model_modifications, analysis)
    for mod in model_modifications
        model = mod(model)
    end
    analysis(model)
end
