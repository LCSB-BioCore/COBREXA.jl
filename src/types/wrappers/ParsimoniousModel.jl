
"""
$(TYPEDEF)

A wrapper that adds a quadratic objective that minimizes the sum of a set of
squared variables. If the bounds of the growth rate of the model are fixed, this
corresponds to finding a parsimonious solution (i.e. pFBA). Note, the sense of
the optimization problem should be `MAX_SENSE`, because internally the objective
is negated.

This is used to implement [`parsimonious_flux_balance_analysis`](@ref).

# Example
```
res = model |>
    with_changed_bound("biomass", lower_bound = 0.1) |>
    with_parsimonious_solution(:enzymes) |>
    flux_balance_analysis(Clarabel.Optimizer)
```
"""
struct ParsimoniousModel <: AbstractModelWrapper
    inner::AbstractMetabolicModel
    var_ids::Vector{String}
end

function ParsimoniousModel(model::AbstractMetabolicModel, semantics::Vector{Symbol})
    var_ids = vcat([first(Accessors.Internal.semantics(sem))(model) for sem in semantics]...)
    ParsimoniousModel(model, var_ids)
end

ParsimoniousModel(model::AbstractMetabolicModel, semantic::Symbol) =
    ParsimoniousModel(model, [semantic])

Accessors.unwrap_model(m::ParsimoniousModel) = m.inner

"""
$(TYPEDSIGNATURES)

Return a negative, uniformly weighted, quadratic-only objective representing the
squared sum of `model.var_ids`.
"""
function Accessors.objective(model::ParsimoniousModel)::SparseMat
    obj = spzeros(n_variables(model), n_variables(model) + 1) # + 1 for QP solver formulation

    idxs = indexin(model.var_ids, variables(model))
    j = findfirst(isnothing, idxs)
    isnothing(j) || throw(DomainError(model.var_ids[j], "No variable with this ID."))

    for i in idxs
        obj[i, i] = -1.0 # negative because objective will be maximized
    end
    obj
end

# need to manually inherit these
Accessors.enzyme_variables(model::ParsimoniousModel) = enzyme_variables(model.inner)
Accessors.enzyme_group_variables(model::ParsimoniousModel) =
    enzyme_group_variables(model.inner)
