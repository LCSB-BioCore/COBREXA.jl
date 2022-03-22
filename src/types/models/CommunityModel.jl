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
mutable struct CommunityModel{M} <: MetabolicModel where {M<:MetabolicModel}
    metabolicModel::M
    exchange_rxn_mets::Dict{String,String}
    biomass_rxn::Maybe{String}
    model_names::Array{String}

    CommunityModel(
        metabolicModel::M = StandardModel();
        exchange_rxn_mets = Dict{String,String}(),
        biomass_rxn = nothing,
        model_names = String[],
    ) where {M<:MetabolicModel} =
        new{M}(metabolicModel, exchange_rxn_mets, biomass_rxn, model_names)
end

@_inherit_model_methods CommunityModel () metabolicModel () reactions metabolites stoichiometry bounds balance objective
