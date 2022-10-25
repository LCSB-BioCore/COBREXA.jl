
"""
$(TYPEDSIGNATURES)

Internal helper to check the presence and shape of modification and argument
arrays in [`screen`](@ref) and pals.
"""
function _screen_args(argtuple, kwargtuple, modsname)

    mods = get(kwargtuple, modsname, nothing)
    args = get(kwargtuple, :args, nothing)

    if isnothing(mods)
        if isnothing(args)
            throw(
                DomainError(args, "at least one of `$modsname` and `args` must be defined"),
            )
        end
        return NamedTuple{(modsname,)}(Ref([[] for _ in args]))
    elseif isnothing(args)
        return (args = [() for _ in mods],)
    else
        if size(mods) != size(args)
            throw(
                DomainError(
                    "$(size(mods)) != $(size(args))",
                    "sizes of `$modsname` and `args` differ",
                ),
            )
        end
        return ()
    end
end

"""
$(TYPEDSIGNATURES)

Take an array of model-modifying function vectors in `variants`, and execute
the function `analysis` on all variants of the `model` specified by `variants`.
The computation is distributed over worker IDs in `workers`. If `args` are
supplied (as an array of the same size as the `variants`), they are forwarded
as arguments to the corresponding analysis function calls.

The array of variants must contain vectors of single-parameter functions, these
are applied to model in order. The functions must *not* modify the model, but
rather return a modified copy. The copy should be made as shallow as possible,
to increase memory efficiency of the process. Variant generators that modify
the argument model in-place will cause unpredictable results. Refer to the
definition of [`screen_variant`](@ref) for details.

The function `analysis` will receive a single argument (the modified model),
together with arguments from `args` expanded by `...`. Supply an array of
tuples or vectors to pass in multiple arguments to each function. If the
argument values should be left intact (not expanded to multiple arguments),
they must be wrapped in single-item tuple or vector.

The modification and analysis functions are transferred to `workers` as-is; all
packages required to run them (e.g. the optimization solvers) must be loaded
there. Typically, you want to use the macro `@everywhere using
MyFavoriteSolver` from `Distributed` package for loading the solvers.

# Return value

The results of running `analysis` are collected in to the resulting array, in a
way that preserves the shape of the `variants`, similarly as with `pmap`.

The results of `analysis` function must be serializable, preferably made only
from pure Julia structures, because they may be transferred over the network
between the computation nodes. For that reason, functions that return whole
JuMP models that contain pointers to allocated C structures (such as
[`flux_balance_analysis`](@ref) used with `GLPK` or `Gurobi` otimizers) will
generally not work in this context.

Note: this function is a thin argument-handling wrapper around
[`_screen_impl`](@ref).

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

screen(m,
    variants = [
        [reverse_reaction(5)],
        [reverse_reaction(3), reverse_reaction(6)]
    ],
    analysis = mod -> mod.S[:,3])  # observe the changes in S

screen(m,
    variants = [
        [reverse_reaction(5)],
        [reverse_reaction(3), reverse_reaction(6)]
    ],
    analysis = mod -> flux_balance_analysis_vec(mod, GLPK.Optimizer))
```
"""
screen(args...; kwargs...) =
    _screen_impl(args...; kwargs..., _screen_args(args, kwargs, :variants)...)

"""
$(TYPEDSIGNATURES)

The actual implementation of [`screen`](@ref).
"""
function _screen_impl(
    model::AbstractMetabolicModel;
    variants::Array{V,N},
    analysis,
    args::Array{A,N},
    workers = [myid()],
)::Array where {V<:AbstractVector,A,N}

    asyncmap(fetch, save_at.(workers, :cobrexa_screen_variants_model, Ref(model)))
    asyncmap(fetch, save_at.(workers, :cobrexa_screen_variants_analysis_fn, Ref(analysis)))
    asyncmap(fetch, get_from.(workers, Ref(:($(precache!)(cobrexa_screen_variants_model)))))

    res = pmap(
        (vars, args)::Tuple -> screen_variant(
            (@remote cobrexa_screen_variants_model),
            vars,
            (@remote cobrexa_screen_variants_analysis_fn),
            args,
        ),
        CachingPool(workers),
        zip(variants, args),
    )

    asyncmap(fetch, remove_from.(workers, :cobrexa_screen_variants_model))
    asyncmap(fetch, remove_from.(workers, :cobrexa_screen_variants_analysis_fn))

    return res
end

"""
$(TYPEDSIGNATURES)

Helper function for [`screen`](@ref) that applies all single-argument
functions in `variant` to the `model` (in order from "first" to
"last"), and executes `analysis` on the result.

Can be used to test model variants locally.
"""
function screen_variant(model::AbstractMetabolicModel, variant::Vector, analysis, args = ())
    for fn in variant
        model = fn(model)
    end
    analysis(model, args...)
end

"""
$(TYPEDSIGNATURES)

A shortcut for [`screen`](@ref) that only works with model variants.
"""
screen_variants(model, variants, analysis; workers = [myid()]) =
    screen(model; variants = variants, analysis = analysis, workers = workers)

"""
$(TYPEDSIGNATURES)

A variant of [`optimize_objective`](@ref) directly usable in
[`screen_optmodel_modifications`](@ref).
"""
screen_optimize_objective(_, optmodel)::Maybe{Float64} = optimize_objective(optmodel)

"""
$(TYPEDSIGNATURES)

Internal helper for [`screen_optmodel_modifications`](@ref) that creates the
model and applies the modifications.
"""
function _screen_optmodel_prepare(model, optimizer, common_modifications)
    precache!(model)
    optmodel = make_optimization_model(model, optimizer)
    for mod in common_modifications
        mod(model, optmodel)
    end
    return (model, optmodel)
end

"""
$(TYPEDSIGNATURES)

Internal helper for [`screen_optmodel_modifications`](@ref) that computes one
item of the screening task.
"""
function _screen_optmodel_item((mods, args))
    (model, optmodel) = @remote cobrexa_screen_optmodel_modifications_data
    for mod in mods
        mod(model, optmodel)
    end
    (@remote cobrexa_screen_optmodel_modifications_fn)(model, optmodel, args...)
end

"""
$(TYPEDSIGNATURES)

Screen multiple modifications of the same optimization model.

This function is potentially more efficient than [`screen`](@ref) because it
avoids making variants of the model structure and remaking of the optimization
model. On the other hand, modification functions need to keep the optimization
model in a recoverable state (one that leaves the model usable for the next
modification), which limits the possible spectrum of modifications applied.

Internally, `model` is distributed to `workers` and transformed into the
optimization model using [`make_optimization_model`](@ref).
`common_modifications` are applied to the models at that point. Next, vectors
of functions in `modifications` are consecutively applied, and the result of
`analysis` function called on model are collected to an array of the same
extent as `modifications`. Calls of `analysis` are optionally supplied with
extra arguments from `args` expanded with `...`, just like in [`screen`](@ref).

Both the modification functions (in vectors) and the analysis function here
have 2 base parameters (as opposed to 1 with [`screen`](@ref)): first is the
`model` (carried through as-is), second is the prepared JuMP optimization model
that may be modified and acted upon. As an example, you can use modification
[`change_constraint`](@ref) and analysis [`screen_optimize_objective`](@ref).

Note: This function is a thin argument-handling wrapper around
[`_screen_optmodel_modifications_impl`](@ref).
"""
screen_optmodel_modifications(args...; kwargs...) = _screen_optmodel_modifications_impl(
    args...;
    kwargs...,
    _screen_args(args, kwargs, :modifications)...,
)

"""
$(TYPEDSIGNATURES)

The actual implementation of [`screen_optmodel_modifications`](@ref).
"""
function _screen_optmodel_modifications_impl(
    model::AbstractMetabolicModel,
    optimizer;
    common_modifications::VF = [],
    modifications::Array{V,N},
    args::Array{A,N},
    analysis::Function,
    workers = [myid()],
)::Array where {V<:AbstractVector,VF<:AbstractVector,A,N}

    asyncmap(
        fetch,
        save_at.(
            workers,
            :cobrexa_screen_optmodel_modifications_data,
            Ref(:($_screen_optmodel_prepare($model, $optimizer, $common_modifications))),
        ),
    )
    asyncmap(
        fetch,
        save_at.(workers, :cobrexa_screen_optmodel_modifications_fn, Ref(analysis)),
    )

    res = pmap(_screen_optmodel_item, CachingPool(workers), zip(modifications, args))

    asyncmap(fetch, remove_from.(workers, :cobrexa_screen_optmodel_modifications_data))
    asyncmap(fetch, remove_from.(workers, :cobrexa_screen_optmodel_modifications_fn))

    return res
end
