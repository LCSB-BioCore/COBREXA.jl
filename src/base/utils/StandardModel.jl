
"""
    Base.copy(m::StandardModel)

Shallow copy of a [`StandardModel`](@ref)
"""
Base.copy(m::StandardModel) = StandardModel(m.id, m.reactions, m.metabolites, m.genes)

"""
    Base.copy(r::Reaction)

Shallow copy of a [`Reaction`](@ref)
"""
Base.copy(r::Reaction) = Reaction(
    r.id;
    name = r.name,
    metabolites = r.metabolites,
    lb = r.lb,
    ub = r.ub,
    grr = r.grr,
    subsystem = r.subsystem,
    notes = r.notes,
    annotations = r.annotations,
    objective_coefficient = r.objective_coefficient,
)

"""
    Base.copy(m::Metabolite)

Shallow copy of a [`Metabolite`](@ref)
"""
Base.copy(m::Metabolite) = Metabolite(
    m.id;
    name = m.name,
    formula = m.formula,
    charge = m.charge,
    compartment = m.compartment,
    notes = m.notes,
    annotations = m.annotations,
)

"""
    Base.copy(g::Gene)

Shallow copy of a [`Gene`](@ref)
"""
Base.copy(g::Gene) = Gene(g.id; name = g.name, notes = g.notes, annotations = g.annotations)

"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model 
given `flux_dict` (the solution of a constraint based analysis) of reactions in `model`.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::StandardModel)
    atom_flux = Dict{String,Float64}()
    for (rxn_id, flux) in flux_dict
        if is_boundary(model.reactions[rxn_id])
            for (met, stoich) in model.reactions[rxn_id].metabolites
                adict = get_atoms(model.metabolites[met])
                for (atom, stoich) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux * stoich
                end
            end
        end
    end
    return atom_flux
end

"""
    metabolite_fluxes(flux_dict::Dict{String, Float64}, model::StandardModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or 
produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(flux_dict::Dict{String,Float64}, model::StandardModel)
    S = stoichiometry(model)
    met_flux = Dict{String,Float64}()
    rxnids = reactions(model)
    metids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = flux_dict[rxnid] * S[row, col]
            # ignore zero flux
            if mf < -_constants.tolerance # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > _constants.tolerance
                if haskey(producing, metid)
                    producing[metid][rxnid] = mf
                else
                    producing[metid] = Dict(rxnid => mf)
                end
            end
        end
    end
    return consuming, producing
end

"""
    set_bound(index, optimization_model;
        ub=_constants.default_reaction_rate,
        lb=-_constants.default_reaction_rate)

Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing, 
so this function simplifies setting constraints. In short, JuMP
uses a normalized right hand side representation of constraints, 
which means that lower bounds have their sign flipped. This function
does this for you, so you don't have to remember to do this whenever you
change the constraints. 

Just supply the constraint `index` and the JuMP model (`opt_model`) that 
will be solved, and the variable's bounds will be set to `ub` and `lb`.
"""
function set_bound(
    vind,
    opt_model;
    ub = _constants.default_reaction_rate,
    lb = -_constants.default_reaction_rate,
)
    set_normalized_rhs(opt_model[:lbs][vind], -lb)
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end

"""
    function find_exchange_reactions(
        model::StandardModel;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )::Vector{String}

Return a vector of reaction ids of exchange reactions. Identifies these reactions based
on their prefixes, see `_constants.exchange_prefixes` for the default list (includes 
prefixes like "EX_", etc.). If `exclude_biomass` is true then the biomass reaction is
not returned as in this list, otherwise looks for the biomass reaction by checking if 
a reaction contains a string in the list `biomass_strings`, see `_constants.biomass_strings`
for the default list (includes strings like "biomass", etc.).

Note: biomass exchange reactions are counted as exchange reactions and will NOT
be excluded when `exclude_biomass` is true.
"""
function find_exchange_reactions(
    model::StandardModel;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Vector{String}
    ex_rxn_ids = String[]
    for rxn_id in reactions(model)
        if any([startswith(rxn_id, x) for x in ex_prefixes]) # found exchange reaction
            push!(ex_rxn_ids, rxn_id)
            continue
        elseif !exclude_biomass && any([occursin(x, rxn_id) for x in biomass_strings]) # biomass
            push!(ex_rxn_ids, rxn_id)    
        end
    end

    return ex_rxn_ids
end

"""
    find_exchange_metabolites(
        model::StandardModel;
        exclude_biomass = false,
        biomass_strings =  _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )::Dict{String, Dict{String, Float64}}

Return a dictionary mapping exchange reaction ids to exchange metabolites in a
dictionary. Identifies these reactions based on their prefixes, see
`_constants.exchange_prefixes` for the default list (includes prefixes like
"EX_", etc.). If `exclude_biomass` is true then the biomass reaction is not
returned as in this list, otherwise looks for the biomass reaction by checking
if a reaction contains a string in the list `biomass_strings`, see
`_constants.biomass_strings` for the default list (includes strings like
"biomass", etc.).

Note: biomass exchange reactions are counted as exchange reactions and will NOT
be excluded when `exclude_biomass` is true.
"""
function find_exchange_metabolites(
    model::StandardModel;
    exclude_biomass = false,
    biomass_strings =  _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Dict{String, Dict{String, Float64}}
    ex_rxns = find_exchange_reactions(
        model,
        exclude_biomass = exclude_biomass,
        biomass_strings = biomass_strings,
        ex_prefixes = ex_prefixes,
    )
    ex_rxn_met = Dict{String, Dict{String, Float64}}()
    for ex_rxn in ex_rxns
        ex_rxn_met[ex_rxn] = model.reactions[ex_rxn].metabolites
    end
    return ex_rxn_met
end
