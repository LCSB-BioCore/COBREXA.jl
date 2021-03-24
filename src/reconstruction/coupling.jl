"""
Add constraints of the following form to a CoupledLinearModel and return a modified one.

The arguments are same as for in-place `addCouplingConstraints!`.
"""
function addCouplingConstraints(m::CoupledLinearModel, args...)
    newLp = deepcopy(m)
    addCouplingConstraints!(newLp, args...)
    return newLp
end

"""
Add constraints to a plain `LinearModel` (converts it to the coupled type)
"""
addCouplingConstraints(m::LinearModel, args...) =
    addCouplingConstraints(convert(CoupledLinearModel, m), args...)

"""
In-place add coupling constraints in form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
function addCouplingConstraints!(
    m::CoupledLinearModel,
    c::V,
    cl::AbstractFloat,
    cu::AbstractFloat,
) where {V<:VecType}
    return addCouplingConstraints!(
        m,
        sparse(reshape(c, (1, length(c)))),
        sparse([cl]),
        sparse([cu]),
    )
end


function addCouplingConstraints!(
    m::CoupledLinearModel,
    C::M,
    cl::V,
    cu::V,
) where {M<:MatType,V<:VecType}

    all([length(cu), length(cl)] .== size(C, 1)) ||
        throw(DimensionMismatch("mismatched numbers of constraints"))
    size(C, 2) == nReactions(m) ||
        throw(DimensionMismatch("mismatched number of reactions"))

    m.C = vcat(m.C, sparse(C))
    m.cl = vcat(m.cl, sparse(cl))
    m.cu = vcat(m.cu, sparse(cu))
end


"""
Remove coupling constraints from the linear model and return the modified model.

Arguments are the same as for in-place version `removeCouplingConstraints!`.
"""
function removeCouplingConstraints(m::CoupledLinearModel, args...)
    newModel = deepcopy(m)
    removeCouplingConstraints!(newModel, args...)
    return newModel
end


"""
Removes a set of coupling constraints from a CoupledLinearModel in-place.
"""
function removeCouplingConstraints!(m::CoupledLinearModel, constraint::Int)
    removeCouplingConstraints!(m, [constraint])
end


function removeCouplingConstraints!(m::CoupledLinearModel, constraints::Vector{Int})
    toBeKept = filter(e -> e ∉ constraints, 1:nCouplingConstraints(m))
    m.C = m.C[toBeKept, :]
    m.cl = m.cl[toBeKept]
    m.cu = m.cu[toBeKept]
end


"""
Change the lower and/or upper bounds ('cl' and 'cu') for given coupling constraints
"""
function changeCouplingBounds!(
    model::CoupledLinearModel,
    constraints::Vector{Int};
    cl::V = Array{Float64}(undef, 0),
    cu::V = Array{Float64}(undef, 0),
) where {V<:VecType}
    found = [index ∈ 1:nCouplingConstraints(model) for index in constraints]
    redConstraints = constraints[found]

    length(redConstraints) == length(unique(redConstraints)) ||
        error("`constraints` appears to contain duplicates")
    if !isempty(cl)
        length(constraints) == length(cl) ||
            throw(DimensionMismatch("`constraints` size doesn't match with `cl`"))
        model.cl[redConstraints] = cl[found]
    end

    if !isempty(cu)
        length(constraints) == length(cu) ||
            throw(DimensionMismatch("`constraints` size doesn't match with `cu`"))
        model.cu[redConstraints] = cu[found]
    end
end
