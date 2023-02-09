
"""
$(TYPEDSIGNATURES)

Pipe-able version of the [`MinimizeAdjustmentModel`](@ref) wrapper.
"""
minimize_adjustment(reference_flux::Vector{Float64}) =
    model -> MinimizeAdjustmentModel(model, reference_flux)
