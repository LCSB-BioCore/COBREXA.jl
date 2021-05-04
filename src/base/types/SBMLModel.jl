"""
    struct SBMLModel

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.
"""
struct SBMLModel <: MetabolicModel
    m::SBML.Model
end

reactions(a::SBMLModel)::Vector{String} = [k for k in keys(a.m.reactions)]
metabolites(a::SBMLModel)::Vector{String} = [k for k in keys(a.m.species)]
n_reactions(a::SBMLModel)::Int = length(a.m.reactions)
n_metabolites(a::SBMLModel)::Int = length(a.m.species)

"""
    stoichiometry(a::SBMLModel)::SparseMat

Recreate the stoichiometry matrix from the SBML model.
"""
function stoichiometry(a::SBMLModel)::SparseMat
    _, _, S = SBML.getS(a.m)
    return S
end

"""
    bounds(a::SBMLModel)::Tuple{SparseVec,SparseVec}

Get the lower and upper flux bounds of a `SBMLModel`. Throws `DomainError` in
case if the SBML contains mismatching units.
"""
function bounds(a::SBMLModel)::Tuple{SparseVec,SparseVec}
    lbu = SBML.getLBs(a.m)
    ubu = SBML.getUBs(a.m)

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

balance(a::SBMLModel)::SparseVec = spzeros(n_metabolites(a))
objective(a::SBMLModel)::SparseVec = SBML.getOCs(a.m)

genes(a::SBMLModel)::Vector{String} = [k for k in a.m.gene_products]
n_genes(a::SBMLModel)::Int = length(a.m.gene_products)

reaction_gene_association(a::SBMLModel, rid::String)::Maybe{GeneAssociation} =
    _maybemap(_parse_grr, a.m.reactions[rid].gene_product_association)

metabolite_formula(a::SBMLModel, mid::String)::Maybe{MetaboliteFormula} =
    _maybemap(_formula_to_atoms, a.m.species[mid].formula)

metabolite_charge(a::SBMLModel, mid::String)::Maybe{Int} = a.m.species[mid].charge

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

_sbml_export_notes = _sbml_export_annotation

function Base.convert(::Type{SBMLModel}, m::MetabolicModel)
    if typeof(m) == SBMLModel
        return m
    end

    mets = metabolites(m)
    rxns = reactions(m)
    stoi = stoichiometry(m)
    (lbs, ubs) = bounds(m)
    ocs = objective(m)
    comps = _default.("", metabolite_compartment.(Ref(m), mets))
    compss = Set(comps)

    return SBMLModel(
        SBML.Model(
            Dict(), # parameters
            Dict("" => []), # units
            Dict(
                [comp =>
                    SBML.Compartment(nothing, nothing, nothing, nothing, nothing, nothing) for
                 comp in compss],
            ),
            Dict(
                [mid => SBML.Species(
                    nothing, # name
                    _default("", comps[mi]), # compartment
                    nothing, # no information about boundary conditions
                    metabolite_formula(m, mid),
                    metabolite_charge(m, mid),
                    nothing, # initial amount
                    nothing, # only substance unit flags
                    _sbml_export_notes(metabolite_notes(m, mid)),
                    _sbml_export_annotation(metabolite_annotations(m, mid)),
                ) for (mi, mid) in enumerate(mets)],
            ),
            Dict(
                [rid => SBML.Reaction(
                    Dict(
                        [mets[i] => stoi[i, ri] for
                         i in SparseArrays.nonzeroinds(stoi[:, ri])],
                    ),
                    (lbs[ri], ""),
                    (ubs[ri], ""),
                    ocs[ri],
                    _maybemap(
                        x -> _unparse_grr(SBML.GeneProductAssociation, x),
                        reaction_gene_association(m, rid),
                    ),
                    nothing,
                    _sbml_export_notes(reaction_notes(m, rid)),
                    _sbml_export_annotation(reaction_annotations(m, rid)),
                ) for (ri, rid) in enumerate(rxns)],
            ),
            Dict(
                [gid => SBML.GeneProduct(
                    nothing,
                    nothing,
                    _sbml_export_notes(gene_notes(m, gid)),
                    _sbml_export_annotation(gene_annotations(m, gid)),
                ) for gid in genes(m)],
            ),
            Dict(), # function definitions
        ),
    )
end
