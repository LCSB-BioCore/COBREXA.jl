
"""
$(TYPEDEF)

A wrapper that adds a quadratic objective that minimizes total squared error
("adjustment") from a certain reference variable assignment ("reference flux").

Length of `reference_assignment` must be the same as the the number of
variables in the inner model.

This is used to implement [`minimize_metabolic_adjustment`](@ref).

# Example
```
m = load_model("e_coli_core.xml")
adjusted_model = m |> MinimizeAdjustmentModel(fill(0.1, n_variables(m)))
```
"""
struct MinimizeAdjustmentModel <: AbstractModelWrapper
    inner::AbstractMetabolicModel
    reference_assignment::Vector{Float64}
end
#TODO: sparse MinimizeAdjustmentModel with dict/sparseVec
#TODO: MinimizeReactionAdjustmentModel?
#TODO: MinimizeEnzymeAdjustmentModel?

Accessors.unwrap_model(m::MinimizeAdjustmentModel) = m.inner

Accessors.objective(m::MinimizeAdjustmentModel)::SparseMat =
    COBREXA.Utils.negative_squared_distance_objective(m.reference_assignment)
