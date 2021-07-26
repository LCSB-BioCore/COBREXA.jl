"""
    struct SBMLModel

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.
"""
struct SBMLModel <: MetabolicModel
    sbml::SBML.Model
end

"""
    reactions(model::SBMLModel)::Vector{String}

Get reactions from a [`SBMLModel`](@ref).
"""
reactions(model::SBMLModel)::Vector{String} = [k for k in keys(model.sbml.reactions)]

"""
    metabolites(model::SBMLModel)::Vector{String}

Get metabolites from a [`SBMLModel`](@ref).
"""
metabolites(model::SBMLModel)::Vector{String} = [k for k in keys(model.sbml.species)]

"""
    n_reactions(model::SBMLModel)::Int

Efficient counting of reactions in [`SBMLModel`](@ref).
"""
n_reactions(model::SBMLModel)::Int = length(model.sbml.reactions)

"""
    n_metabolites(model::SBMLModel)::Int

Efficient counting of metabolites in [`SBMLModel`](@ref).
"""
n_metabolites(model::SBMLModel)::Int = length(model.sbml.species)

"""
    stoichiometry(model::SBMLModel)::SparseMat

Recreate the stoichiometry matrix from the [`SBMLModel`](@ref).
"""
function stoichiometry(model::SBMLModel)::SparseMat
    _, _, S = SBML.getS(model.sbml)
    return S
end

"""
    bounds(model::SBMLModel)::Tuple{SparseVec,SparseVec}

Get the lower and upper flux bounds of model [`SBMLModel`](@ref). Throws `DomainError` in
case if the SBML contains mismatching units.
"""
function bounds(model::SBMLModel)::Tuple{SparseVec,SparseVec}
    lbu = SBML.getLBs(model.sbml)
    ubu = SBML.getUBs(model.sbml)

    unit = lbu[1][2]
    getvalue = (val, _)::Tuple -> val
    getunit = (_, unit)::Tuple -> unit

    allunits = unique([getunit.(lbu) getunit.(ubu)])
    length(allunits) == 1 || throw(
        DomainError(
            allunits,
            "The SBML file uses multiple units; loading needs conversion",
        ),
    )

    return sparse.((getvalue.(lbu), getvalue.(ubu)))
end

"""
    balance(model::SBMLModel)::SparseVec

Balance vector of a [`SBMLModel`](@ref). This is always zero.
"""
balance(model::SBMLModel)::SparseVec = spzeros(n_metabolites(model))

"""
    objective(model::SBMLModel)::SparseVec

Objective of the [`SBMLModel`](@ref).
"""
objective(model::SBMLModel)::SparseVec = SBML.getOCs(model.sbml)

"""
    genes(model::SBMLModel)::Vector{String}

Get genes of a [`SBMLModel`](@ref).
"""
genes(model::SBMLModel)::Vector{String} = [k for k in keys(model.sbml.gene_products)]

"""
    n_genes(model::SBMLModel)::Int

Get number of genes in [`SBMLModel`](@ref).
"""
n_genes(model::SBMLModel)::Int = length(model.sbml.gene_products)

"""
    reaction_gene_association(model::SBMLModel, rid::String)::Maybe{GeneAssociation}

Retrieve the [`GeneAssociation`](@ref) from [`SBMLModel`](@ref).
"""
reaction_gene_association(model::SBMLModel, rid::String)::Maybe{GeneAssociation} =
    _maybemap(_parse_grr, model.sbml.reactions[rid].gene_product_association)

"""
    metabolite_formula(model::SBMLModel, mid::String)::Maybe{MetaboliteFormula}

Get [`MetaboliteFormula`](@ref) from a chosen metabolite from [`SBMLModel`](@ref).
"""
metabolite_formula(model::SBMLModel, mid::String)::Maybe{MetaboliteFormula} =
    _maybemap(_parse_formula, model.sbml.species[mid].formula)

"""
    metabolite_charge(model::SBMLModel, mid::String)::Maybe{Int}

Get charge of a chosen metabolite from [`SBMLModel`](@ref).
"""
metabolite_charge(model::SBMLModel, mid::String)::Maybe{Int} =
    model.sbml.species[mid].charge

function _sbml_export_annotation(annotation)::Maybe{String}
    if isnothing(annotation) || isempty(annotation)
        nothing
    elseif length(annotation) != 1 || first(annotation).first != ""
        @_io_log @warn "Possible data loss: annotation cannot be exported to SBML" annotation
        join(["$k: $v" for (k, v) in annotation], "\n")
    else
        first(annotation).second
    end
end

const _sbml_export_notes = _sbml_export_annotation

"""
    reaction_stoichiometry(model::SBMLModel, rid::String)::Dict{String, Float64}

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::SBMLModel, rid::String)::Dict{String,Float64} =
    m.sbml.reactions[rid].stoichiometry

"""
    Base.convert(::Type{SBMLModel}, mm::MetabolicModel)

Convert any metabolic model to [`SBMLModel`](@ref).
"""
function Base.convert(::Type{SBMLModel}, mm::MetabolicModel)
    if typeof(mm) == SBMLModel
        return mm
    end

    mets = metabolites(mm)
    rxns = reactions(mm)
    stoi = stoichiometry(mm)
    (lbs, ubs) = bounds(mm)
    ocs = objective(mm)
    comps = _default.("", metabolite_compartment.(Ref(mm), mets))
    compss = Set(comps)

    return SBMLModel(
        SBML.Model(
            Dict(), # parameters
            Dict("" => []), # units
            Dict([
                comp =>
                    SBML.Compartment(nothing, nothing, nothing, nothing, nothing, nothing)
                for comp in compss
            ],),
            Dict([
                mid => SBML.Species(
                    nothing, # name
                    _default("", comps[mi]), # compartment
                    nothing, # no information about boundary conditions
                    metabolite_formula(mm, mid),
                    metabolite_charge(mm, mid),
                    nothing, # initial amount
                    nothing, # initial concentration
                    nothing, # only substance unit flags
                    _sbml_export_notes(metabolite_notes(mm, mid)),
                    _sbml_export_annotation(metabolite_annotations(mm, mid)),
                ) for (mi, mid) in enumerate(mets)
            ],),
            Dict([
                rid => SBML.Reaction(
                    Dict([
                        mets[i] => stoi[i, ri] for
                        i in SparseArrays.nonzeroinds(stoi[:, ri])
                    ],),
                    (lbs[ri], ""),
                    (ubs[ri], ""),
                    ocs[ri],
                    _maybemap(
                        x -> _unparse_grr(SBML.GeneProductAssociation, x),
                        reaction_gene_association(mm, rid),
                    ),
                    nothing, # no kinetic math
                    true, # reversible by default
                    _sbml_export_notes(reaction_notes(mm, rid)),
                    _sbml_export_annotation(reaction_annotations(mm, rid)),
                ) for (ri, rid) in enumerate(rxns)
            ],),
            Dict([
                gid => SBML.GeneProduct(
                    nothing,
                    nothing,
                    _sbml_export_notes(gene_notes(mm, gid)),
                    _sbml_export_annotation(gene_annotations(mm, gid)),
                ) for gid in genes(mm)
            ],),
            Dict(), # function definitions
        ),
    )
end
