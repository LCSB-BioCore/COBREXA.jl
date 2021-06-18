"""
    add_reactions(
        m::CoreModelCoupled,
        s::V1,
        b::V2,
        c::AbstractFloat,
        xl::AbstractFloat,
        xu::AbstractFloat;
        check_consistency = false,
    ) where {V1<:VecType,V2<:VecType}

Add reaction(s) to a `CoreModelCoupled` model `m`.

"""
function add_reactions(
    m::CoreModelCoupled,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType}
    new_lm = add_reactions(m.lm, s, b, c, xl, xu, check_consistency = check_consistency)
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm) - n_reactions(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
    add_reactions(
        m::CoreModelCoupled,
        s::V1,
        b::V2,
        c::AbstractFloat,
        xl::AbstractFloat,
        xu::AbstractFloat,
        rxn::String,
        mets::K;
        check_consistency = false,
    ) where {V1<:VecType,V2<:VecType,K<:StringVecType}

"""
function add_reactions(
    m::CoreModelCoupled,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat,
    rxn::String,
    mets::K;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType,K<:StringVecType}
    new_lm = add_reactions(
        m.lm,
        s,
        b,
        c,
        xl,
        xu,
        rxn,
        mets,
        check_consistency = check_consistency,
    )
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm) - n_reactions(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
    add_reactions(
        m::CoreModelCoupled,
        Sp::M,
        b::V,
        c::V,
        xl::V,
        xu::V;
        check_consistency = false,
    ) where {M<:MatType,V<:VecType}

"""
function add_reactions(
    m::CoreModelCoupled,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V;
    check_consistency = false,
) where {M<:MatType,V<:VecType}
    new_lm = add_reactions(m.lm, Sp, b, c, xl, xu, check_consistency = check_consistency)
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm) - n_reactions(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
    add_reactions(m1::CoreModelCoupled, m2::CoreModel; check_consistency = false)

Add all reactions from `m2` to `m1`.

"""
function add_reactions(m1::CoreModelCoupled, m2::CoreModel; check_consistency = false)
    new_lm = add_reactions(m1.lm, m2, check_consistency = check_consistency)
    return CoreModelCoupled(
        new_lm,
        hcat(m1.C, spzeros(size(m1.C, 1), n_reactions(new_lm) - n_reactions(m1.lm))),
        m1.cl,
        m1.cu,
    )
end

"""
    add_reactions(
        m::CoreModelCoupled,
        Sp::M,
        b::V,
        c::V,
        xl::V,
        xu::V,
        rxns::K,
        mets::K;
        check_consistency = false,
    ) where {M<:MatType,V<:VecType,K<:StringVecType}

"""
function add_reactions(
    m::CoreModelCoupled,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    rxns::K,
    mets::K;
    check_consistency = false,
) where {M<:MatType,V<:VecType,K<:StringVecType}
    new_lm = add_reactions(
        m.lm,
        Sp,
        b,
        c,
        xl,
        xu,
        rxns,
        mets,
        check_consistency = check_consistency,
    )
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm) - n_reactions(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
    remove_reactions(m::CoreModelCoupled, rxns::Vector{Int})

Remove reaction(s) from a `CoreModelCoupled`.

Also removes any metabolites not involved in any reaction after the deletion.

"""
function remove_reactions(m::CoreModelCoupled, rxns::Vector{Int})
    return CoreModelCoupled(
        remove_reactions(m.lm, rxns),
        m.C[:, filter(e -> e ∉ rxns, 1:n_reactions(m))],
        m.cl,
        m.cu,
    )
end

"""
    remove_reactions(m::CoreModelCoupled, rxn::Integer)

"""
function remove_reactions(m::CoreModelCoupled, rxn::Integer)
    return remove_reactions(m, [rxn])
end

"""
    remove_reactions(m::CoreModelCoupled, rxn::String)

"""
function remove_reactions(m::CoreModelCoupled, rxn::String)
    return remove_reactions(m, [rxn])
end

"""
    remove_reactions(m::CoreModelCoupled, rxns::Vector{String})

"""
function remove_reactions(m::CoreModelCoupled, rxns::Vector{String})
    rxn_indices =
        [findfirst(isequal(name), m.lm.rxns) for name in intersect(rxns, m.lm.rxns)]
    if isempty(rxn_indices)
        return m
    else
        return remove_reactions(m, rxn_indices)
    end
end

"""
Add constraints of the following form to a CoreModelCoupled and return a modified one.

Add constraints to a [`CoreModelCoupled`](@ref) and return a modified one.

The arguments are same as for in-place [`add_coupling_constraints!`](@ref).
"""
function add_coupling_constraints(m::CoreModelCoupled, args...)
    new_lp = deepcopy(m)
    add_coupling_constraints!(new_lp, args...)
    return new_lp
end

"""
    add_coupling_constraints(m::CoreModel, args...)

Add coupling constraints to a plain [`CoreModel`](@ref) (returns a
[`CoreModelCoupled`](@ref)).
"""
add_coupling_constraints(m::CoreModel, args...) =
    add_coupling_constraints(convert(CoreModelCoupled, m), args...)

"""
    add_coupling_constraints!(
        m::CoreModelCoupled,
        c::VecType,
        cl::AbstractFloat,
        cu::AbstractFloat,
    )

Overload for adding a single coupling constraint.
"""
function add_coupling_constraints!(
    m::CoreModelCoupled,
    c::VecType,
    cl::AbstractFloat,
    cu::AbstractFloat,
)
    return add_coupling_constraints!(
        m,
        sparse(reshape(c, (1, length(c)))),
        sparse([cl]),
        sparse([cu]),
    )
end


"""
    add_coupling_constraints!(
        m::CoreModelCoupled,
        C::MatType,
        cl::V,
        cu::V,
    ) where {V<:VecType}

In-place add a single coupling constraint in form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
function add_coupling_constraints!(
    m::CoreModelCoupled,
    C::MatType,
    cl::V,
    cu::V,
) where {V<:VecType}

    all([length(cu), length(cl)] .== size(C, 1)) ||
        throw(DimensionMismatch("mismatched numbers of constraints"))
    size(C, 2) == n_reactions(m) ||
        throw(DimensionMismatch("mismatched number of reactions"))

    m.C = vcat(m.C, sparse(C))
    m.cl = vcat(m.cl, sparse(cl))
    m.cu = vcat(m.cu, sparse(cu))
end


"""
    remove_coupling_constraints(m::CoreModelCoupled, args...)

Remove coupling constraints from the linear model, and return the modified
model. Arguments are the same as for in-place version
[`remove_coupling_constraints!`](@ref).
"""
function remove_coupling_constraints(m::CoreModelCoupled, args...)
    new_model = deepcopy(m)
    remove_coupling_constraints!(new_model, args...)
    return new_model
end


"""
    remove_coupling_constraints!(m::CoreModelCoupled, constraint::Int)

Removes a single coupling constraints from a [`CoreModelCoupled`](@ref)
in-place.
"""
remove_coupling_constraints!(m::CoreModelCoupled, constraint::Int) =
    remove_coupling_constraints!(m, [constraint])


"""
    remove_coupling_constraints!(m::CoreModelCoupled, constraints::Vector{Int})

Removes a set of coupling constraints from a [`CoreModelCoupled`](@ref)
in-place.
"""
function remove_coupling_constraints!(m::CoreModelCoupled, constraints::Vector{Int})
    to_be_kept = filter(e -> e ∉ constraints, 1:n_coupling_constraints(m))
    m.C = m.C[to_be_kept, :]
    m.cl = m.cl[to_be_kept]
    m.cu = m.cu[to_be_kept]
end


"""
    change_coupling_bounds!(
        model::CoreModelCoupled,
        constraints::Vector{Int};
        cl::V = Float64[],
        cu::V = Float64[],
    ) where {V<:VecType}

Change the lower and/or upper bounds (`cl` and `cu`) for the given list of
coupling constraints.
"""
function change_coupling_bounds!(
    model::CoreModelCoupled,
    constraints::Vector{Int};
    cl::V = Float64[],
    cu::V = Float64[],
) where {V<:VecType}
    found = [index ∈ 1:n_coupling_constraints(model) for index in constraints]
    red_constraints = constraints[found]

    length(red_constraints) == length(unique(red_constraints)) ||
        error("`constraints` appears to contain duplicates")
    if !isempty(cl)
        length(constraints) == length(cl) ||
            throw(DimensionMismatch("`constraints` size doesn't match with `cl`"))
        model.cl[red_constraints] = cl[found]
    end

    if !isempty(cu)
        length(constraints) == length(cu) ||
            throw(DimensionMismatch("`constraints` size doesn't match with `cu`"))
        model.cu[red_constraints] = cu[found]
    end
end

"""
    change_bounds!(
        model::CoreModelCoupled,
        reaction_idxs::Vector{Int};
        lower_bounds = fill(-_constants.default_reaction_bound, length(reaction_idxs)),
        upper_bounds = fill(constants.default_reaction_bound, length(reaction_idxs)),
    )   

Change the bounds of all reactions with indices `reaction_idxs` in `model` in-place. Note
that if the bound argument is not supplied then a default (see
`_constants.default_reaction_bound`) is used. The same default bound argument is
used for each reaction bound if not supplied. 

See also: [`change_bound`](@ref), [`change_bound!`](@ref), [`change_bounds`](@ref) 

# Example
```
change_bounds!(model, [2, 3]; lower_bounds=[-10, -20], ub=[10, 22])
change_bounds!(model, [2, 3]; lower_bounds=[-10.2, -14.3]) # all upper_bounds are set to _constants.default_reaction_bound
change_bounds!(model, [2, 3]; upper_bounds=[10.2, 23])
```
"""
function change_bounds!(
    model::CoreModelCoupled,
    reaction_idxs::Vector{Int};
    lower_bounds = fill(-_constants.default_reaction_bound, length(rxns)),
    upper_bounds = fill(_constants.default_reaction_bound, length(rxns)),
)   
    for (rxn_idx, lb, ub) in zip(reaction_idxs, lower_bounds, upper_bounds)
        change_bound!(model, rxn_idx; lower_bound=lb, upper_bound=ub)
    end
end

@_change_bound!("CoreModelCoupled",)
function change_bound!(
    model::CoreModelCoupled,
    rxn::Int;
    lower_bound = -_constants.default_reaction_bound,
    upper_bound = _constants.default_reaction_bound,
)
    model.lm.xl[rxn] = lower_bound    
    model.lm.xu[rxn] = upper_bound
    return nothing # so that nothing gets printed
end

"""
    change_bounds!(
        model::CoreModelCoupled,
        rxn_ids::Vector{String};
        lower_bounds = fill(-_constants.default_reaction_bound, length(rxn_ids)),
        upper_bounds = fill(constants.default_reaction_bound, length(rxn_ids)),
    )

Change the bounds of all reactions with `rxn_ids` in `model` in-place. Note
that if the bound argument is not supplied then a default (see
`_constants.default_reaction_bound`) is used. The same default bound argument is
used for each reaction bound if not supplied. 

See also: [`change_bound`](@ref), [`change_bound!`](@ref), [`change_bounds`](@ref) 

# Example
```
change_bounds!(model, ["PFL", "FBA"]; lower_bounds=[-10, -20], ub=[10, 22])
change_bounds!(model, ["PFL", "FBA"]; lower_bounds=[-10.2, -14.3]) # all upper_bounds are set to _constants.default_reaction_bound
change_bounds!(model, ["PFL", "FBA"]; upper_bounds=[10.2, 23])
```
"""
function change_bounds!(
    model::CoreModelCoupled,
    rxn_ids::Vector{String};
    lower_bounds = fill(-_constants.default_reaction_bound, length(rxn_ids)),
    upper_bounds = fill(constants.default_reaction_bound, length(rxn_ids)),
)
    change_bounds!(model, Int.(indexin(rxn_ids, reactions(model))); lower_bounds = lower_bounds, upper_bounds = upper_bounds)
end

@_change_bound!("CoreModelCoupled",)
function change_bound!(
    model::CoreModelCoupled,
    rxn_id::String;
    lower_bound = -_constants.default_reaction_bound,
    upper_bound = _constants.default_reaction_bound,
)
    change_bound!(model, first(indexin([rxn_id], reactions(model))); lower_bound = lower_bound, upper_bound = upper_bound)
end

"""
    change_bounds(
        model::CoreModelCoupled,
        rxns::Vector{Int};
        lower_bounds = fill(-_constants.default_reaction_bound, length(rxns)),
        upper_bounds = fill(constants.default_reaction_bound, length(rxns)),
    ) 

Return a shallow copy of the `model` where the bounds of all reactions with
indices `rxns` in `model` have been changed. Note that if the bound
argument is not supplied then a default (see
`_constants.default_reaction_bound`) is used. The same default bound argument is
used for each reaction bound if not supplied. 

See also: [`change_bound`](@ref), [`change_bound!`](@ref), [`change_bounds`](@ref) 

# Example
```
new_model = change_bounds(model, [2, 3]; lower_bounds=[-10, -20], ub=[10, 22])
new_model = change_bounds(model, [2, 3]; lower_bounds=[-10.2, -14.3]) # all upper_bounds are set to _constants.default_reaction_bound
new_model = change_bounds(model, [2, 3]; upper_bounds=[10.2, 23])
```
"""
function change_bounds(
    model::CoreModelCoupled,
    rxns::Vector{Int};
    lower_bounds = fill(-_constants.default_reaction_bound, length(rxns)),
    upper_bounds = fill(constants.default_reaction_bound, length(rxns)),
)       
    m = copy(model)
    m.lm.xl = copy(model.lm.xl)
    m.lm.xu = copy(model.lm.xu)
    m.lm.xl[rxns] .= lower_bounds    
    m.lm.xu[rxns] .= upper_bounds
    return m
end

"""
    change_bound(
        model::CoreModelCoupled,
        reaction_idx::Int;
        lower_bound = -_constants.default_reaction_bound,
        upper_bound = _constants.default_reaction_bound,
    )

Return a shallow copy of the `model` where the bound of a reaction with index
`reaction_idx` in `model` has been changed. Note that if the bound argument is
not supplied then a default (see `_constants.default_reaction_bound`) is used.
The same default bound argument is used for each reaction bound if not supplied.

See also: [`change_bound`](@ref), [`change_bound!`](@ref), [`change_bounds`](@ref) 

# Example
```
new_model = change_bound(model, 2; lower_bounds=-10, ub=22)
new_model = change_bound(model, 2; lower_bounds=-10.2) # all upper_bounds are set to _constants.default_reaction_bound
new_model = change_bound(model, 2; upper_bounds=10.2)
```
"""
function change_bound(
    model::CoreModelCoupled,
    reaction_idx::Int;
    lower_bound = -_constants.default_reaction_bound,
    upper_bound = _constants.default_reaction_bound,
)
    m = copy(model)
    m.lm.xl = copy(model.lm.xl)
    m.lm.xu = copy(model.lm.xu)
    m.lm.xl[reaction_idx] = lower_bound
    m.lm.xu[reaction_idx] = upper_bound
    return m
end

"""
    change_bounds(
        model::CoreModelCoupled,
        rxn_ids::Vector{String};
        lower_bounds = fill(-_constants.default_reaction_bound, length(rxn_ids)),
        upper_bounds = fill(constants.default_reaction_bound, length(rxn_ids)),
    )

Change the bounds of all reactions with `rxn_ids` in `model` in-place. Note
that if the bound argument is not supplied then a default (see
`_constants.default_reaction_bound`) is used. The same default bound argument is
used for each reaction bound if not supplied. 

See also: [`change_bound`](@ref), [`change_bound!`](@ref), [`change_bounds`](@ref) 

# Example
```
new_model = change_bounds(model, ["PFL", "FBA"]; lower_bounds=[-10, -20], ub=[10, 22])
new_model = change_bounds(model, ["PFL", "FBA"]; lower_bounds=[-10.2, -14.3]) # all upper_bounds are set to _constants.default_reaction_bound
new_model = change_bounds(model, ["PFL", "FBA"]; upper_bounds=[10.2, 23])
```
"""
function change_bounds(
    model::CoreModelCoupled,
    rxn_ids::Vector{String};
    lower_bounds = fill(-_constants.default_reaction_bound, length(rxn_ids)),
    upper_bounds = fill(constants.default_reaction_bound, length(rxn_ids)),
)
    change_bounds(model, Int.(indexin(rxn_ids, reactions(model))); lower_bounds = lower_bounds, upper_bounds = upper_bounds)
end

"""
    change_bound(
        model::CoreModelCoupled,
        reaction_id::String;
        lower_bound = -_constants.default_reaction_bound,
        upper_bound = _constants.default_reaction_bound,
    )

Change the bounds of a reaction with `reaction_id` in `model` in-place. Note
that if the bound argument is not supplied then a default (see
`_constants.default_reaction_bound`) is used. 

See also: [`change_bound`](@ref), [`change_bounds!`](@ref), [`change_bounds!`](@ref) 

# Example
```
new_model = change_bound!(model, "PFL"; lower_bound=-10, ub=10)
new_model = change_bound!(model, "PFL"; lower_bound=-10.2) # upper_bound is set to _constants.default_reaction_bound
new_model = change_bound!(model, "PFL"; upper_bound=10)
```
"""
function change_bound(
    model::CoreModelCoupled,
    rxn_id::String;
    lower_bound = -_constants.default_reaction_bound,
    upper_bound = _constants.default_reaction_bound,
)
    change_bound(model, first(indexin([rxn_id], reactions(model))); lower_bound = lower_bound, upper_bound = upper_bound)
end
