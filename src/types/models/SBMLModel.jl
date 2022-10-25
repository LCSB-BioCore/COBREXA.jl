"""
$(TYPEDEF)

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.

# Fields
$(TYPEDFIELDS)
"""
struct SBMLModel <: AbstractMetabolicModel
    sbml::SBML.Model
end

"""
$(TYPEDSIGNATURES)

Get reactions from a [`SBMLModel`](@ref).
"""
Accessors.reactions(model::SBMLModel)::Vector{String} =
    [k for k in keys(model.sbml.reactions)]

"""
$(TYPEDSIGNATURES)

Get metabolites from a [`SBMLModel`](@ref).
"""
Accessors.metabolites(model::SBMLModel)::Vector{String} =
    [k for k in keys(model.sbml.species)]

"""
$(TYPEDSIGNATURES)

Efficient counting of reactions in [`SBMLModel`](@ref).
"""
Accessors.n_reactions(model::SBMLModel)::Int = length(model.sbml.reactions)

"""
$(TYPEDSIGNATURES)

Efficient counting of metabolites in [`SBMLModel`](@ref).
"""
Accessors.n_metabolites(model::SBMLModel)::Int = length(model.sbml.species)

"""
$(TYPEDSIGNATURES)

Recreate the stoichiometry matrix from the [`SBMLModel`](@ref).
"""
function Accessors.stoichiometry(model::SBMLModel)::SparseMat
    _, _, S = SBML.stoichiometry_matrix(model.sbml)
    return S
end

"""
$(TYPEDSIGNATURES)

Get the lower and upper flux bounds of model [`SBMLModel`](@ref). Throws `DomainError` in
case if the SBML contains mismatching units.
"""
function Accessors.bounds(model::SBMLModel)::Tuple{Vector{Float64},Vector{Float64}}
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
Accessors.balance(model::SBMLModel)::SparseVec = spzeros(n_metabolites(model))

"""
$(TYPEDSIGNATURES)

Objective of the [`SBMLModel`](@ref).
"""
Accessors.objective(model::SBMLModel)::SparseVec = SBML.flux_objective(model.sbml)

"""
$(TYPEDSIGNATURES)

Get genes of a [`SBMLModel`](@ref).
"""
Accessors.genes(model::SBMLModel)::Vector{String} =
    [k for k in keys(model.sbml.gene_products)]

"""
$(TYPEDSIGNATURES)

Get number of genes in [`SBMLModel`](@ref).
"""
Accessors.n_genes(model::SBMLModel)::Int = length(model.sbml.gene_products)

"""
$(TYPEDSIGNATURES)

Retrieve the [`GeneAssociation`](@ref) from [`SBMLModel`](@ref).
"""
Accessors.reaction_gene_association(model::SBMLModel, rid::String)::Maybe{GeneAssociation} =
    maybemap(parse_grr, model.sbml.reactions[rid].gene_product_association)

"""
$(TYPEDSIGNATURES)

Get [`MetaboliteFormula`](@ref) from a chosen metabolite from [`SBMLModel`](@ref).
"""
Accessors.metabolite_formula(model::SBMLModel, mid::String)::Maybe{MetaboliteFormula} =
    maybemap(parse_formula, model.sbml.species[mid].formula)

"""
$(TYPEDSIGNATURES)

Get the compartment of a chosen metabolite from [`SBMLModel`](@ref).
"""
Accessors.metabolite_compartment(model::SBMLModel, mid::String) =
    model.sbml.species[mid].compartment

"""
$(TYPEDSIGNATURES)

Get charge of a chosen metabolite from [`SBMLModel`](@ref).
"""
Accessors.metabolite_charge(model::SBMLModel, mid::String)::Maybe{Int} =
    model.sbml.species[mid].charge

function _parse_sbml_identifiers_org_uri(uri::String)::Tuple{String,String}
    m = match(r"^http://identifiers.org/([^/]+)/(.*)$", uri)
    isnothing(m) ? ("RESOURCE_URI", uri) : (m[1], m[2])
end

function _sbml_import_cvterms(sbo::Maybe{String}, cvs::Vector{SBML.CVTerm})::Annotations
    res = Annotations()
    isnothing(sbo) || (res["sbo"] = [sbo])
    for cv in cvs
        cv.biological_qualifier == :is || continue
        for (id, val) in _parse_sbml_identifiers_org_uri.(cv.resource_uris)
            push!(get!(res, id, []), val)
        end
    end
    return res
end

function _sbml_export_cvterms(annotations::Annotations)::Vector{SBML.CVTerm}
    isempty(annotations) && return []
    length(annotations) == 1 && haskey(annotations, "sbo") && return []
    [
        SBML.CVTerm(
            biological_qualifier = :is,
            resource_uris = [
                id == "RESOURCE_URI" ? val : "http://identifiers.org/$id/$val" for
                (id, vals) = annotations if id != "sbo" for val in vals
            ],
        ),
    ]
end

function _sbml_export_sbo(annotations::Annotations)::Maybe{String}
    haskey(annotations, "sbo") || return nothing
    if length(annotations["sbo"]) != 1
        @io_log @error "Data loss: SBO term is not unique for SBML export" annotations["sbo"]
        return
    end
    return annotations["sbo"][1]
end

function _sbml_import_notes(notes::Maybe{String})::Notes
    isnothing(notes) ? Notes() : Notes("" => [notes])
end

function _sbml_export_notes(notes::Notes)::Maybe{String}
    isempty(notes) || @io_log @error "Data loss: notes not exported to SBML" notes
    nothing
end

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
"""
function Accessors.reaction_stoichiometry(m::SBMLModel, rid::String)::Dict{String,Float64}
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
Accessors.reaction_name(model::SBMLModel, rid::String) = model.sbml.reactions[rid].name

"""
$(TYPEDSIGNATURES)

Return the name of metabolite with ID `mid`.
"""
Accessors.metabolite_name(model::SBMLModel, mid::String) = model.sbml.species[mid].name

"""
$(TYPEDSIGNATURES)

Return the name of gene with ID `gid`.
"""
Accessors.gene_name(model::SBMLModel, gid::String) = model.sbml.gene_products[gid].name

"""
$(TYPEDSIGNATURES)

Return the annotations of reaction with ID `rid`.
"""
Accessors.reaction_annotations(model::SBMLModel, rid::String) =
    _sbml_import_cvterms(model.sbml.reactions[rid].sbo, model.sbml.reactions[rid].cv_terms)

"""
$(TYPEDSIGNATURES)

Return the annotations of metabolite with ID `mid`.
"""
Accessors.metabolite_annotations(model::SBMLModel, mid::String) =
    _sbml_import_cvterms(model.sbml.species[mid].sbo, model.sbml.species[mid].cv_terms)

"""
$(TYPEDSIGNATURES)

Return the annotations of gene with ID `gid`.
"""
Accessors.gene_annotations(model::SBMLModel, gid::String) = _sbml_import_cvterms(
    model.sbml.gene_products[gid].sbo,
    model.sbml.gene_products[gid].cv_terms,
)

"""
$(TYPEDSIGNATURES)

Return the notes about reaction with ID `rid`.
"""
Accessors.reaction_notes(model::SBMLModel, rid::String) =
    _sbml_import_notes(model.sbml.reactions[rid].notes)

"""
$(TYPEDSIGNATURES)

Return the notes about metabolite with ID `mid`.
"""
Accessors.metabolite_notes(model::SBMLModel, mid::String) =
    _sbml_import_notes(model.sbml.species[mid].notes)

"""
$(TYPEDSIGNATURES)

Return the notes about gene with ID `gid`.
"""
Accessors.gene_notes(model::SBMLModel, gid::String) =
    _sbml_import_notes(model.sbml.gene_products[gid].notes)

"""
$(TYPEDSIGNATURES)

Convert any metabolic model to [`SBMLModel`](@ref).
"""
function Base.convert(::Type{SBMLModel}, mm::AbstractMetabolicModel)
    if typeof(mm) == SBMLModel
        return mm
    end

    mets = metabolites(mm)
    rxns = reactions(mm)
    stoi = stoichiometry(mm)
    (lbs, ubs) = bounds(mm)
    comps = default.("compartment", metabolite_compartment.(Ref(mm), mets))
    compss = Set(comps)

    metid(x) = startswith(x, "M_") ? x : "M_$x"
    rxnid(x) = startswith(x, "R_") ? x : "R_$x"
    gprid(x) = startswith(x, "G_") ? x : "G_$x"

    return SBMLModel(
        SBML.Model(
            compartments = Dict(
                comp => SBML.Compartment(constant = true) for comp in compss
            ),
            species = Dict(
                metid(mid) => SBML.Species(
                    name = metabolite_name(mm, mid),
                    compartment = default("compartment", comps[mi]),
                    formula = maybemap(unparse_formula, metabolite_formula(mm, mid)),
                    charge = metabolite_charge(mm, mid),
                    constant = false,
                    boundary_condition = false,
                    only_substance_units = false,
                    sbo = _sbml_export_sbo(metabolite_annotations(mm, mid)),
                    notes = _sbml_export_notes(metabolite_notes(mm, mid)),
                    metaid = metid(mid),
                    cv_terms = _sbml_export_cvterms(metabolite_annotations(mm, mid)),
                ) for (mi, mid) in enumerate(mets)
            ),
            reactions = Dict(
                rxnid(rid) => SBML.Reaction(
                    name = reaction_name(mm, rid),
                    reactants = [
                        SBML.SpeciesReference(
                            species = metid(mets[i]),
                            stoichiometry = -stoi[i, ri],
                            constant = true,
                        ) for
                        i in SparseArrays.nonzeroinds(stoi[:, ri]) if stoi[i, ri] <= 0
                    ],
                    products = [
                        SBML.SpeciesReference(
                            species = metid(mets[i]),
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
                    gene_product_association = maybemap(
                        x -> unparse_grr(SBML.GeneProductAssociation, x),
                        reaction_gene_association(mm, rid),
                    ),
                    reversible = true,
                    sbo = _sbml_export_sbo(reaction_annotations(mm, rid)),
                    notes = _sbml_export_notes(reaction_notes(mm, rid)),
                    metaid = rxnid(rid),
                    cv_terms = _sbml_export_cvterms(reaction_annotations(mm, rid)),
                ) for (ri, rid) in enumerate(rxns)
            ),
            gene_products = Dict(
                gprid(gid) => SBML.GeneProduct(
                    label = gid,
                    name = gene_name(mm, gid),
                    sbo = _sbml_export_sbo(gene_annotations(mm, gid)),
                    notes = _sbml_export_notes(gene_notes(mm, gid)),
                    metaid = gprid(gid),
                    cv_terms = _sbml_export_cvterms(gene_annotations(mm, gid)),
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
