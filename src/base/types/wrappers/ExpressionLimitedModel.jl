
"""
    $(TYPEDEF)

[`ExpressionLimitedModel`](@ref) follows the methodology of the E-Flux algorithm to
constraint the flux through the reactions in order to simulate the limited
expression of genes (and thus the limited availability of the gene products).

Use [`make_expression_limited_model`](@ref) or [`with_expression_limits`](@ref)
to construct the models.

E-Flux algorithm is closer described by: *Colijn, Caroline, Aaron Brandes,
Jeremy Zucker, Desmond S. Lun, Brian Weiner, Maha R. Farhat, Tan-Yun Cheng, D.
Branch Moody, Megan Murray, and James E. Galagan. "Interpreting expression data
with metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid
production." PLoS computational biology 5, no. 8 (2009): e1000489*.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct ExpressionLimitedModel <: ModelWrapper
    """
    Relative gene expression w.r.t. to some chosen reference; the normalization
    and scale of the values should match the expectations of the
    `bounding_function`.
    """
    relative_expression::Dict{String,Float64}

    "The wrapped model."
    inner::MetabolicModel

    """
    Function used to calculate the new reaction bounds from expression.
    In [`make_expression_limited_model`](@ref) this defaults to
    [`expression_probabilistic_bounds`](@ref).
    """
    bounding_function::Function
end

COBREXA.unwrap_model(m::ExpressionLimitedModel) = m.inner

function COBREXA.bounds(m::ExpressionLimitedModel)::Tuple{Vector{Float64},Vector{Float64}}
    (lbs, ubs) = bounds(m.inner)
    lims = collect(
        m.bounding_function(m.relative_expression, reaction_gene_association(m.inner, rid)) for rid in reactions(m.inner)
    )
    (lbs .* lims, ubs .* lims)
end
