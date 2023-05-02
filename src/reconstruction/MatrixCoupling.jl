"""
$(TYPEDSIGNATURES)

Add reaction(s) to a `MatrixModelWithCoupling` model `m`.

"""
function add_reactions(
    m::MatrixModelWithCoupling,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType}
    new_lm = add_reactions(m.lm, s, b, c, xl, xu, check_consistency = check_consistency)
    return MatrixModelWithCoupling(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), variable_count(new_lm) - variable_count(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
$(TYPEDSIGNATURES)

"""
function add_reactions(
    m::MatrixModelWithCoupling,
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
    return MatrixModelWithCoupling(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), variable_count(new_lm) - variable_count(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
$(TYPEDSIGNATURES)
"""
function add_reactions(
    m::MatrixModelWithCoupling,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V;
    check_consistency = false,
) where {M<:MatType,V<:VecType}
    new_lm = add_reactions(m.lm, Sp, b, c, xl, xu, check_consistency = check_consistency)
    return MatrixModelWithCoupling(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), variable_count(new_lm) - variable_count(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
$(TYPEDSIGNATURES)

Add all reactions from `m2` to `m1`.
"""
function add_reactions(
    m1::MatrixModelWithCoupling,
    m2::MatrixModel;
    check_consistency = false,
)
    new_lm = add_reactions(m1.lm, m2, check_consistency = check_consistency)
    return MatrixModelWithCoupling(
        new_lm,
        hcat(m1.C, spzeros(size(m1.C, 1), variable_count(new_lm) - variable_count(m1.lm))),
        m1.cl,
        m1.cu,
    )
end

"""
$(TYPEDSIGNATURES)
"""
function add_reactions(
    m::MatrixModelWithCoupling,
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
    return MatrixModelWithCoupling(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), variable_count(new_lm) - variable_count(m.lm))),
        m.cl,
        m.cu,
    )
end

"""
$(TYPEDSIGNATURES)

Add constraints of the following form to MatrixCoupling and return the modified
model.

The arguments are same as for in-place [`add_coupling_constraints!`](@ref).
"""
function add_coupling_constraints(m::MatrixCoupling, args...)
    new_lp = deepcopy(m)
    add_coupling_constraints!(new_lp, args...)
    return new_lp
end

"""
$(TYPEDSIGNATURES)

Add coupling constraints to a plain [`MatrixModel`](@ref) (returns a
[`MatrixModelWithCoupling`](@ref)).
"""
add_coupling_constraints(m::MatrixModel, args...) = MatrixModelWithCoupling(m, args...)

"""
$(TYPEDSIGNATURES)

Overload for adding a single coupling constraint.
"""
function add_coupling_constraints!(
    m::MatrixCoupling,
    c::VecType,
    cl::AbstractFloat,
    cu::AbstractFloat,
)
    return add_coupling_constraints!(m, sparse(reshape(c, (1, length(c)))), [cl], [cu])
end

"""
$(TYPEDSIGNATURES)

In-place add a single coupling constraint in form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
function add_coupling_constraints!(
    m::MatrixCoupling,
    C::MatType,
    cl::V,
    cu::V,
) where {V<:VecType}

    all([length(cu), length(cl)] .== size(C, 1)) ||
        throw(DimensionMismatch("mismatched numbers of constraints"))
    size(C, 2) == variable_count(m) ||
        throw(DimensionMismatch("mismatched number of reactions"))

    m.C = vcat(m.C, sparse(C))
    m.cl = vcat(m.cl, collect(cl))
    m.cu = vcat(m.cu, collect(cu))
    nothing
end

"""
$(TYPEDSIGNATURES)

Remove coupling constraints from the linear model, and return the modified
model. Arguments are the same as for in-place version
[`remove_coupling_constraints!`](@ref).
"""
function remove_coupling_constraints(m::MatrixCoupling, args...)
    new_model = deepcopy(m)
    remove_coupling_constraints!(new_model, args...)
    return new_model
end

"""
$(TYPEDSIGNATURES)

Removes a single coupling constraints from a [`MatrixCoupling`](@ref) in-place.
"""
remove_coupling_constraints!(m::MatrixCoupling, constraint::Int) =
    remove_coupling_constraints!(m, [constraint])


"""
$(TYPEDSIGNATURES)

Removes a set of coupling constraints from a [`MatrixCoupling`](@ref)
in-place.
"""
function remove_coupling_constraints!(m::MatrixCoupling, constraints::Vector{Int})
    to_be_kept = filter(!in(constraints), 1:n_coupling_constraints(m))
    m.C = m.C[to_be_kept, :]
    m.cl = m.cl[to_be_kept]
    m.cu = m.cu[to_be_kept]
    nothing
end

"""
$(TYPEDSIGNATURES)

Change the lower and/or upper bounds (`cl` and `cu`) for the given list of
coupling constraints.
"""
function change_coupling_bounds!(
    model::MatrixCoupling,
    constraints::Vector{Int};
    cl::V = Float64[],
    cu::V = Float64[],
) where {V<:VecType}
    found = (constraints .>= 1) .& (constraints .<= n_coupling_constraints(model))
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
    nothing
end

# TODO see if some of these can be derived from AbstractModelWrapper
@_change_bounds_fn MatrixCoupling Int inplace begin
    change_bound!(model.lm, rxn_idx; lower_bound, upper_bound)
end

@_change_bounds_fn MatrixCoupling Int inplace plural begin
    change_bounds!(model.lm, rxn_idxs; lower_bounds, upper_bounds)
end

@_change_bounds_fn MatrixCoupling String inplace begin
    change_bound!(model.lm, rxn_id; lower_bound, upper_bound)
end

@_change_bounds_fn MatrixCoupling String inplace plural begin
    change_bounds!(model.lm, rxn_ids; lower_bounds, upper_bounds)
end

@_change_bounds_fn MatrixCoupling Int begin
    n = copy(model)
    n.lm = change_bound(model.lm, rxn_idx; lower_bound, upper_bound)
    n
end

@_change_bounds_fn MatrixCoupling Int plural begin
    n = copy(model)
    n.lm = change_bounds(model.lm, rxn_idxs; lower_bounds, upper_bounds)
    n
end

@_change_bounds_fn MatrixCoupling String begin
    n = copy(model)
    n.lm = change_bound(model.lm, rxn_id; lower_bound, upper_bound)
    n
end

@_change_bounds_fn MatrixCoupling String plural begin
    n = copy(model)
    n.lm = change_bounds(model.lm, rxn_ids; lower_bounds, upper_bounds)
    n
end

@_remove_fn reaction MatrixCoupling Int inplace begin
    remove_reactions!(model, [reaction_idx])
end

@_remove_fn reaction MatrixCoupling Int inplace plural begin
    orig_rxns = variable_ids(model.lm)
    remove_reactions!(model.lm, reaction_idxs)
    model.C = model.C[:, in.(orig_rxns, Ref(Set(variable_ids(model.lm))))]
    nothing
end

@_remove_fn reaction MatrixCoupling Int begin
    remove_reactions(model, [reaction_idx])
end

@_remove_fn reaction MatrixCoupling Int plural begin
    n = copy(model)
    n.lm = remove_reactions(n.lm, reaction_idxs)
    n.C = n.C[:, in.(variable_ids(model.lm), Ref(Set(variable_ids(n.lm))))]
    return n
end

@_remove_fn reaction MatrixCoupling String inplace begin
    remove_reactions!(model, [reaction_id])
end

@_remove_fn reaction MatrixCoupling String inplace plural begin
    remove_reactions!(model, Int.(indexin(reaction_ids, variable_ids(model))))
end

@_remove_fn reaction MatrixCoupling String begin
    remove_reactions(model, [reaction_id])
end

@_remove_fn reaction MatrixCoupling String plural begin
    remove_reactions(model, Int.(indexin(reaction_ids, variable_ids(model))))
end

@_remove_fn metabolite MatrixCoupling Int inplace begin
    remove_metabolites!(model, [metabolite_idx])
end

@_remove_fn metabolite MatrixCoupling Int plural inplace begin
    orig_rxns = variable_ids(model.lm)
    model.lm = remove_metabolites(model.lm, metabolite_idxs)
    model.C = model.C[:, in.(orig_rxns, Ref(Set(variable_ids(model.lm))))]
    nothing
end

@_remove_fn metabolite MatrixCoupling Int begin
    remove_metabolites(model, [metabolite_idx])
end

@_remove_fn metabolite MatrixCoupling Int plural begin
    n = copy(model)
    n.lm = remove_metabolites(n.lm, metabolite_idxs)
    return n
end

@_remove_fn metabolite MatrixCoupling String inplace begin
    remove_metabolites!(model, [metabolite_id])
end

@_remove_fn metabolite MatrixCoupling String inplace plural begin
    remove_metabolites!(model, Int.(indexin(metabolite_ids, metabolites(model))))
end

@_remove_fn metabolite MatrixCoupling String begin
    remove_metabolites(model, [metabolite_id])
end

@_remove_fn metabolite MatrixCoupling String plural begin
    remove_metabolites(model, Int.(indexin(metabolite_ids, metabolites(model))))
end

"""
$(TYPEDSIGNATURES)

Forwards arguments to [`change_objective!`](@ref) of the internal model.
"""
function change_objective!(model::MatrixCoupling, args...; kwargs...)
    change_objective!(model.lm, args...; kwargs...)
end
