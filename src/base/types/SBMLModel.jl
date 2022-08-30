"""
    struct SBMLModel

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.
"""
struct SBMLModel <: MetabolicModel
    sbml::SBML.Model
end

"""
$(TYPEDSIGNATURES)

Get reactions from a [`SBMLModel`](@ref).
"""
reactions(model::SBMLModel)::Vector{String} = [k for k in keys(model.sbml.reactions)]

"""
$(TYPEDSIGNATURES)

Get metabolites from a [`SBMLModel`](@ref).
"""
metabolites(model::SBMLModel)::Vector{String} = [k for k in keys(model.sbml.species)]

"""
$(TYPEDSIGNATURES)

Efficient counting of reactions in [`SBMLModel`](@ref).
"""
n_reactions(model::SBMLModel)::Int = length(model.sbml.reactions)

"""
$(TYPEDSIGNATURES)

Efficient counting of metabolites in [`SBMLModel`](@ref).
"""
n_metabolites(model::SBMLModel)::Int = length(model.sbml.species)

"""
$(TYPEDSIGNATURES)

Recreate the stoichiometry matrix from the [`SBMLModel`](@ref).
"""
function stoichiometry(model::SBMLModel)::SparseMat
    _, _, S = SBML.stoichiometry_matrix(model.sbml)
    return S
end

"""
$(TYPEDSIGNATURES)

Get the lower and upper flux bounds of model [`SBMLModel`](@ref). Throws `DomainError` in
case if the SBML contains mismatching units.
"""
function bounds(model::SBMLModel)::Tuple{Vector{Float64},Vector{Float64}}
    lbu, ubu = SBML.flux_bounds(model.sbml)

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

    return (getvalue.(lbu), getvalue.(ubu))
end

"""
$(TYPEDSIGNATURES)

Balance vector of a [`SBMLModel`](@ref). This is always zero.
"""
balance(model::SBMLModel)::SparseVec = spzeros(n_metabolites(model))

"""
$(TYPEDSIGNATURES)

Objective of the [`SBMLModel`](@ref).
"""
objective(model::SBMLModel)::SparseVec = SBML.flux_objective(model.sbml)

"""
$(TYPEDSIGNATURES)

Get genes of a [`SBMLModel`](@ref).
"""
genes(model::SBMLModel)::Vector{String} = [k for k in keys(model.sbml.gene_products)]

"""
$(TYPEDSIGNATURES)

Get number of genes in [`SBMLModel`](@ref).
"""
n_genes(model::SBMLModel)::Int = length(model.sbml.gene_products)

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`SBMLModel`](@ref).
"""
reaction_gene_association(model::SBMLModel, rid::String)::Maybe{GeneAssociation} =
    _maybemap(_parse_grr, model.sbml.reactions[rid].gene_product_association)

"""
$(TYPEDSIGNATURES)

Get [`MetaboliteFormula`](@ref) from a chosen metabolite from [`SBMLModel`](@ref).
"""
metabolite_formula(model::SBMLModel, mid::String)::Maybe{MetaboliteFormula} =
    _maybemap(_parse_formula, model.sbml.species[mid].formula)

"""
$(TYPEDSIGNATURES)

Get charge of a chosen metabolite from [`SBMLModel`](@ref).
"""
metabolite_charge(model::SBMLModel, mid::String)::Maybe{Int} =
    model.sbml.species[mid].charge

function _sbml_export_annotation(annotation)::Maybe{String}
    if isnothing(annotation) || isempty(annotation)
        nothing
    elseif length(annotation) != 1 || first(annotation).first != ""
        @_io_log @warn "Possible data loss: multiple annotations converted to text for SBML" annotation
        join(["$k: $v" for (k, v) in annotation], "\n")
    else
        @_io_log @warn "Possible data loss: trying to represent annotation in SBML is unlikely to work " annotation
        first(annotation).second
    end
end

const _sbml_export_notes = _sbml_export_annotation

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
"""
function reaction_stoichiometry(m::SBMLModel, rid::String)::Dict{String,Float64}
    s = Dict{String,Float64}()
    default1(x) = isnothing(x) ? 1 : x
    for sr in m.sbml.reactions[rid].reactants
        s[sr.species] = get(s, sr.species, 0.0) - default1(sr.stoichiometry)
    end
    for sr in m.sbml.reactions[rid].products
        s[sr.species] = get(s, sr.species, 0.0) + default1(sr.stoichiometry)
    end
    return s
end

"""
$(TYPEDSIGNATURES)

Return the name of reaction with ID `rid`.
"""
reaction_name(model::SBMLModel, rid::String) = model.sbml.reactions[rid].name

"""
$(TYPEDSIGNATURES)

Return the name of metabolite with ID `mid`.
"""
metabolite_name(model::SBMLModel, mid::String) = model.sbml.species[mid].name

"""
$(TYPEDSIGNATURES)

Return the name of gene with ID `gid`.
"""
gene_name(model::SBMLModel, gid::String) = model.sbml.gene_products[gid].name

"""
$(TYPEDSIGNATURES)

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
    comps = _default.("compartment", metabolite_compartment.(Ref(mm), mets))
    compss = Set(comps)

    return SBMLModel(
        SBML.Model(
            compartments = Dict(
                comp => SBML.Compartment(constant = true) for comp in compss
            ),
            species = Dict(
                mid => SBML.Species(
                    name = metabolite_name(mm, mid),
                    compartment = _default("compartment", comps[mi]),
                    formula = metabolite_formula(mm, mid),
                    charge = metabolite_charge(mm, mid),
                    constant = false,
                    boundary_condition = false,
                    only_substance_units = false,
                    notes = _sbml_export_notes(metabolite_notes(mm, mid)),
                    annotation = _sbml_export_annotation(metabolite_annotations(mm, mid)),
                ) for (mi, mid) in enumerate(mets)
            ),
            reactions = Dict(
                rid => SBML.Reaction(
                    name = reaction_name(mm, rid),
                    reactants = [
                        SBML.SpeciesReference(
                            species = mets[i],
                            stoichiometry = -stoi[i, ri],
                            constant = true,
                        ) for
                        i in SparseArrays.nonzeroinds(stoi[:, ri]) if stoi[i, ri] <= 0
                    ],
                    products = [
                        SBML.SpeciesReference(
                            species = mets[i],
                            stoichiometry = stoi[i, ri],
                            constant = true,
                        ) for
                        i in SparseArrays.nonzeroinds(stoi[:, ri]) if stoi[i, ri] > 0
                    ],
                    kinetic_parameters = Dict(
                        "LOWER_BOUND" => SBML.Parameter(value = lbs[ri]),
                        "UPPER_BOUND" => SBML.Parameter(value = ubs[ri]),
                    ),
                    lower_bound = "LOWER_BOUND",
                    upper_bound = "UPPER_BOUND",
                    gene_product_association = _maybemap(
                        x -> _unparse_grr(SBML.GeneProductAssociation, x),
                        reaction_gene_association(mm, rid),
                    ),
                    reversible = true,
                    notes = _sbml_export_notes(reaction_notes(mm, rid)),
                    annotation = _sbml_export_annotation(reaction_annotations(mm, rid)),
                ) for (ri, rid) in enumerate(rxns)
            ),
            gene_products = Dict(
                gid => SBML.GeneProduct(
                    label = gid,
                    name = gene_name(mm, gid),
                    notes = _sbml_export_notes(gene_notes(mm, gid)),
                    annotation = _sbml_export_annotation(gene_annotations(mm, gid)),
                ) for gid in genes(mm)
            ),
            active_objective = "objective",
            objectives = Dict(
                "objective" => SBML.Objective(
                    "maximize",
                    Dict(rid => oc for (rid, oc) in zip(rxns, objective(mm)) if oc != 0),
                ),
            ),
        ),
    )
end
