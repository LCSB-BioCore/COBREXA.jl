
"""
$(TYPEDSIGNATURES)

Get the internal [`CoreModel`](@ref) out of [`CoreCoupling`](@ref).
"""
unwrap_model(a::CoreCoupling) = a.lm

"""
$(TYPEDSIGNATURES)

Coupling constraint matrix for a `CoreCoupling`.
"""
coupling(a::CoreCoupling)::SparseMat = vcat(coupling(a.lm), a.C)

"""
$(TYPEDSIGNATURES)

The number of coupling constraints in a `CoreCoupling`.
"""
n_coupling_constraints(a::CoreCoupling)::Int = n_coupling_constraints(a.lm) + size(a.C, 1)

"""
$(TYPEDSIGNATURES)

Coupling bounds for a `CoreCoupling`.
"""
coupling_bounds(a::CoreCoupling)::Tuple{Vector{Float64},Vector{Float64}} =
    vcat.(coupling_bounds(a.lm), (a.cl, a.cu))

"""
$(TYPEDSIGNATURES)

Make a `CoreCoupling` out of any compatible model type.
"""
function Base.convert(
    ::Type{CoreCoupling{M}},
    mm::MetabolicModel;
    clone_coupling = true,
) where {M}
    if mm isa CoreCoupling{M}
        mm
    elseif mm isa CoreCoupling
        CoreCoupling(convert(M, mm.lm), mm.C, mm.cl, mm.cu)
    elseif clone_coupling
        (cl, cu) = coupling_bounds(mm)
        CoreCoupling(convert(M, mm), coupling(mm), cl, cu)
    else
        CoreCoupling(convert(M, mm), spzeros(0, n_reactions(mm)), spzeros(0), spzeros(0))
    end
end

# these are special for CoreModel-ish models
@_inherit_model_methods CoreModelCoupled () lm () reaction_gene_association_vec
@_inherit_model_methods CoreModelCoupled (ridx::Int,) lm (ridx,) reaction_stoichiometry
