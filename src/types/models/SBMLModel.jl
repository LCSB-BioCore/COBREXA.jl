"""
$(TYPEDEF)

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.

# Fields
$(TYPEDFIELDS)
"""
struct SBMLModel <: AbstractMetabolicModel
    sbml::SBML.Model
    reaction_ids::Vector{String}
    reaction_idx::Dict{String,Int}
    metabolite_ids::Vector{String}
    metabolite_idx::Dict{String,Int}
    gene_ids::Vector{String}
    active_objective::String
end

"""
$(TYPEDEF)

Construct the SBML model and add the necessary cached indexes, possibly choosing an active objective.
"""
function SBMLModel(sbml::SBML.Model, active_objective::String = "")
    rxns = sort(collect(keys(sbml.reactions)))
    mets = sort(collect(keys(sbml.species)))
    genes = sort(collect(keys(sbml.gene_products)))

    SBMLModel(
        sbml,
        rxns,
        Dict(rxns .=> eachindex(rxns)),
        mets,
        Dict(mets .=> eachindex(mets)),
        genes,
        active_objective,
    )
end

Accessors.variable_ids(model::SBMLModel)::Vector{String} = model.reaction_ids

Accessors.Internal.@all_variables_are_reactions SBMLModel

Accessors.metabolite_ids(model::SBMLModel)::Vector{String} = model.metabolite_ids

function Accessors.metabolite_variables(model::SBMLModel)
    x = Dict{String,Dict{String,Float64}}()
    add!(mid, rid, c) =
        if haskey(x, mid)
            x[mid][rid] = c
        else
            x[mid] = Dict{String,Float64}(rid => c)
        end
    for (rid, r) in model.sbml.reactions
        for reactant in r.reactants
            add!(
                reactant.species,
                rid,
                isnothing(reactant.stoichiometry) ? -1.0 : -reactant.stoichiometry,
            )
        end
        for product in r.products
            add!(
                product.species,
                rid,
                isnothing(product.stoichiometry) ? 1.0 : product.stoichiometry,
            )
        end
    end
    x
end

function Accessors.variable_bounds(model::SBMLModel)::Tuple{Vector{Float64},Vector{Float64}}
    # There are multiple ways in SBML to specify a lower/upper bound. There are
    # the "global" model bounds that we completely ignore now because no one
    # uses them. In reaction, you can specify the bounds using "LOWER_BOUND"
    # and "UPPER_BOUND" parameters, but also there may be a FBC plugged-in
    # parameter name that refers to the parameters.  We extract these, using
    # the units from the parameters. For unbounded reactions we use -Inf or Inf
    # as a default.

    common_unit = ""

    function get_bound(rid, fld, param, default)
        rxn = model.sbml.reactions[rid]
        param_name = SBML.mayfirst(getfield(rxn, fld), param)
        param = get(
            rxn.kinetic_parameters,
            param_name,
            get(model.sbml.parameters, param_name, default),
        )
        unit = SBML.mayfirst(param.units, "")
        if unit != ""
            if common_unit != ""
                if unit != common_unit
                    throw(
                        DomainError(
                            units_in_sbml,
                            "The SBML file uses multiple units; loading would need conversion",
                        ),
                    )
                end
            else
                common_unit = unit
            end
        end
        return param.value
    end

    return (
        get_bound.(
            model.reaction_ids,
            :lower_bound,
            "LOWER_BOUND",
            Ref(SBML.Parameter(value = -Inf)),
        ),
        get_bound.(
            model.reaction_ids,
            :upper_bound,
            "UPPER_BOUND",
            Ref(SBML.Parameter(value = Inf)),
        ),
    )
end

Accessors.metabolite_bounds(model::SBMLModel)::SparseVec = spzeros(metabolite_count(model))

function Accessors.objective(model::SBMLModel)::SparseVec
    res = sparsevec([], [], reaction_count(model))

    objective = get(model.sbml.objectives, model.active_objective, nothing)
    if isnothing(objective) && length(model.sbml.objectives) == 1
        objective = first(values(model.sbml.objectives))
    end
    if !isnothing(objective)
        direction = objective.type == "maximize" ? 1.0 : -1.0
        for (rid, coef) in objective.flux_objectives
            res[model.reaction_idx[rid]] = float(direction * coef)
        end
    else
        # old-style objectives
        for (rid, r) in model.sbml.reactions
            oc = get(r.kinetic_parameters, "OBJECTIVE_COEFFICIENT", nothing)
            isnothing(oc) || (res[model.reaction_idx[rid]] = float(oc.value))
        end
    end
    return res
end

Accessors.genes(model::SBMLModel)::Vector{String} = model.gene_ids

Accessors.reaction_gene_associations(
    model::SBMLModel,
    rid::String,
)::Maybe{GeneAssociationsDNF} =
    maybemap(parse_grr, model.sbml.reactions[rid].gene_product_association)

Accessors.eval_reaction_gene_association(model::SBMLModel, rid::String; kwargs...) =
    maybemap(
        x -> eval_grr(x; kwargs...),
        model.sbml.reactions[rid].gene_product_association,
    )

Accessors.metabolite_formula(model::SBMLModel, mid::String)::Maybe{MetaboliteFormula} =
    maybemap(parse_formula, model.sbml.species[mid].formula)

Accessors.metabolite_compartment(model::SBMLModel, mid::String) =
    model.sbml.species[mid].compartment

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
                (id, vals) in annotations if id != "sbo" for val in vals
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

Accessors.reaction_name(model::SBMLModel, rid::String) = model.sbml.reactions[rid].name

Accessors.metabolite_name(model::SBMLModel, mid::String) = model.sbml.species[mid].name

Accessors.gene_name(model::SBMLModel, gid::String) = model.sbml.gene_products[gid].name

Accessors.reaction_annotations(model::SBMLModel, rid::String) =
    _sbml_import_cvterms(model.sbml.reactions[rid].sbo, model.sbml.reactions[rid].cv_terms)

Accessors.metabolite_annotations(model::SBMLModel, mid::String) =
    _sbml_import_cvterms(model.sbml.species[mid].sbo, model.sbml.species[mid].cv_terms)

Accessors.gene_annotations(model::SBMLModel, gid::String) = _sbml_import_cvterms(
    model.sbml.gene_products[gid].sbo,
    model.sbml.gene_products[gid].cv_terms,
)

Accessors.reaction_notes(model::SBMLModel, rid::String) =
    _sbml_import_notes(model.sbml.reactions[rid].notes)

Accessors.metabolite_notes(model::SBMLModel, mid::String) =
    _sbml_import_notes(model.sbml.species[mid].notes)

Accessors.gene_notes(model::SBMLModel, gid::String) =
    _sbml_import_notes(model.sbml.gene_products[gid].notes)

Accessors.model_annotations(model::SBMLModel) =
    _sbml_import_cvterms(model.sbml.sbo, model.sbml.cv_terms)

Accessors.model_notes(model::SBMLModel) = _sbml_import_notes(model.sbml.notes)

function Base.convert(::Type{SBMLModel}, mm::AbstractMetabolicModel)
    if typeof(mm) == SBMLModel
        return mm
    end

    mets = metabolite_ids(mm)
    rxns = variable_ids(mm)
    stoi = stoichiometry(mm)
    (lbs, ubs) = variable_bounds(mm)
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
                        reaction_gene_associations(mm, rid),
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
            notes = _sbml_export_notes(model_notes(mm)),
            sbo = _sbml_export_sbo(model_annotations(mm)),
            cv_terms = _sbml_export_cvterms(model_annotations(mm)),
        ),
    )
end
