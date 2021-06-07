"""
    add_coupling_constraints(m::CoreModelCoupled, args...)

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
