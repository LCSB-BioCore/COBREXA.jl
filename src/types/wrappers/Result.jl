"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
struct Result{ResultType} <: AbstractModelWrapper
    solved_model::ResultType
    model::AbstractMetabolicModel
end

Accessors.unwrap_model(result::Result) = result.model