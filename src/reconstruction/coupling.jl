"""
Add constraints of the following form to a CoupledLinearModel and return a modified one.

The arguments are same as for in-place `add_coupling_constraints!`.
"""
function add_coupling_constraints(m::CoupledLinearModel, args...)
    new_lp = deepcopy(m)
    add_coupling_constraints!(new_lp, args...)
    return new_lp
end

"""
Add constraints to a plain `LinearModel` (converts it to the coupled type)
"""
add_coupling_constraints(m::LinearModel, args...) =
    add_coupling_constraints(convert(CoupledLinearModel, m), args...)

"""
In-place add coupling constraints in form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
function add_coupling_constraints!(
    m::CoupledLinearModel,
    c::V,
    cl::AbstractFloat,
    cu::AbstractFloat,
) where {V<:VecType}
    return add_coupling_constraints!(
        m,
        sparse(reshape(c, (1, length(c)))),
        sparse([cl]),
        sparse([cu]),
    )
end


function add_coupling_constraints!(
    m::CoupledLinearModel,
    C::M,
    cl::V,
    cu::V,
) where {M<:MatType,V<:VecType}

    all([length(cu), length(cl)] .== size(C, 1)) ||
        throw(DimensionMismatch("mismatched numbers of constraints"))
    size(C, 2) == n_reactions(m) ||
        throw(DimensionMismatch("mismatched number of reactions"))

    m.C = vcat(m.C, sparse(C))
    m.cl = vcat(m.cl, sparse(cl))
    m.cu = vcat(m.cu, sparse(cu))
end


"""
Remove coupling constraints from the linear model and return the modified model.

Arguments are the same as for in-place version `remove_coupling_constraints!`.
"""
function remove_coupling_constraints(m::CoupledLinearModel, args...)
    new_model = deepcopy(m)
    remove_coupling_constraints!(new_model, args...)
    return new_model
end


"""
Removes a set of coupling constraints from a CoupledLinearModel in-place.
"""
function remove_coupling_constraints!(m::CoupledLinearModel, constraint::Int)
    remove_coupling_constraints!(m, [constraint])
end


function remove_coupling_constraints!(m::CoupledLinearModel, constraints::Vector{Int})
    to_be_kept = filter(e -> e ∉ constraints, 1:n_coupling_constraints(m))
    m.C = m.C[to_be_kept, :]
    m.cl = m.cl[to_be_kept]
    m.cu = m.cu[to_be_kept]
end


"""
Change the lower and/or upper bounds ('cl' and 'cu') for given coupling constraints
"""
function change_coupling_bounds!(
    model::CoupledLinearModel,
    constraints::Vector{Int};
    cl::V = Array{Float64}(undef, 0),
    cu::V = Array{Float64}(undef, 0),
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
