"""
    mutable struct CommunityModel

`CommunityModel` is used to store a so-called compartmentalised community model
where a number of individual models are joined together via a joined
compartment, which in turn connects the whole model to the "outside".
Additionally, the model may include a community level biomass (pseudo)reaction.

The actual model is stored as a regular `MetabolicModel` that functions just
like a "normal" model would. In addition, `exchange_rxn_mets` stores the
community level exchanges connecting the individual models together
(see [`join_with_exchanges`](@ref)), `biomass_rxn` is the id of the
community biomass reaction (if it exists), and `model_names` lists the string
identifiers of the individual models. It is expected that `model_names`
corresponds to the prefixes used to identify which reaction belongs to which
model.

See also: [`join_with_exchanges`](@ref)

# Fields
```
metabolicModel :: MetabolicModel
exchange_rxn_mets :: Dict{String,String}
biomass_rxn :: Maybe{String}
model_names :: Array{String}
```

# Example
```

```
"""
mutable struct CommunityModel <: MetabolicModel
    metabolicModel::Any
    exchange_rxn_mets::Dict{String,String}
    biomass_rxn::Maybe{String}
    model_names::Array{String}

    CommunityModel(
        metabolicModel = StandardModel();
        exchange_rxn_mets = Dict{String,String}(),
        biomass_rxn = nothing,
        model_names = String[],
    ) = new(metabolicModel, exchange_rxn_mets, biomass_rxn, model_names)
end

"""
    reactions(cm::CommunityModel)

Get the reactions from the `MetabolicModel` in a `CommunityModel`.
"""
reactions(cm::CommunityModel) = reactions(cm.metabolicModel)

"""
    metabolites(cm::CommunityModel)

Get the metabolites from the `MetabolicModel` in a `CommunityModel`.
"""
metabolites(cm::CommunityModel) = metabolites(cm.metabolicModel)

"""
    stoichiometry(cm::CommunityModel)

Get the stoichiometry of the underlying `MetabolicModel`.
"""
stoichiometry(cm::CommunityModel) = stoichiometry(cm.metabolicModel)

"""
    bounds(cm::CommunityModel)

Get the flux bounds of the underlying `MetabolicModel`.
"""
bounds(cm::CommunityModel) = bounds(cm.metabolicModel)

"""
    objective(cm::CommunityModel)

Get the objective vector of the underlying `MetabolicModel`.
"""
objective(cm::CommunityModel) = objective(cm.metabolicModel)
