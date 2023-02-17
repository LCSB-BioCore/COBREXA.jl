"""
    mutable struct Result

A `Result` type stores an [`AbstractMetabolicModel`](@ref) and a `JuMP.Model`. 
"""
mutable struct Result{AM<:AbstractMetabolicModel,JM<:JuMP.Model} <: AbstractResult
    model::AM
    opt_model::JM
end
