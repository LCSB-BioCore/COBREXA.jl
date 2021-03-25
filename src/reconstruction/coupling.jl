"""
Verifies that vectors and matrices have the expected dimensions.
"""
function checkCouplingConstraintsInputDimensions(
    m::LinearModel,
    C::M,
    cl::V,
    cu::V,
) where {M<:MatType,V<:VecType}

    length(cu) == length(cl) ||
        throw(DimensionMismatch("`cl` and `cu` don't have the same size"))
    size(C) == (length(cl), nReactions(m)) || throw(
        DimensionMismatch("`C` size should be length(cl) x the number of reactions in m"),
    )
end

"""
Add constraints of the following form to a LinearModel and return a modified one.

The arguments are same as for in-place `addCouplingConstraints!`.
"""
function addCouplingConstraints(m::LinearModel, args...)
    newLp = deepcopy(m)
    addCouplingConstraints!(newLp, args...)
    return newLp
end

"""
In-place add coupling constraints in form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
function addCouplingConstraints!(
    m::LinearModel,
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
    m::LinearModel,
    C::M,
    cl::V,
    cu::V,
) where {M<:MatType,V<:VecType}

    C = sparse(C)
    cl = sparse(cl)
    cu = sparse(cu)

    checkCouplingConstraintsInputDimensions(m, C, cl, cu)

    m.C = vcat(m.C, C)
    m.cl = vcat(m.cl, cl)
    m.cu = vcat(m.cu, cu)
end


"""
Remove coupling constraints from the linear model and return the modified model.

Arguments are the same as for in-place version `removeCouplingConstraints!`.
"""
function removeCouplingConstraints(m::LinearModel, args...)
    newModel = deepcopy(m)
    removeCouplingConstraints!(newModel, args...)
    return newModel
end


"""
Removes a set of coupling constraints from a LinearModel in-place.
"""
function removeCouplingConstraints!(m::LinearModel, constraint::Int)
    removeCouplingConstraints!(m, [constraint])
end


function removeCouplingConstraints!(m::LinearModel, constraints::Vector{Int})
    toBeKept = filter(e -> e ∉ constraints, 1:nCouplingConstraints(m))
    m.C = m.C[toBeKept, :]
    m.cl = m.cl[toBeKept]
    m.cu = m.cu[toBeKept]
end


"""
Change the lower and/or upper bounds ('cl' and 'cu') for given coupling constraints
"""
function changeCouplingBounds!(
    model::LinearModel,
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
