# # Building and analysing a small community model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we will use `COBREXA` to build and analyze a small community model  
# consisting of multiple variants of *E. coli* using the `CoreModel`. Each
# variant will only be able to metabolize one substrate as its energy/carbon
# source. The environment will supply only sugars, meaning that some of the
# variants will only be able to feed off of the waste products of the sugar
# metabolizing microbes. We will use an objective function that enforces equal
# growth rates.

# ## Load the base model

!isfile("iML1515.json") &&
    download("http://bigg.ucsd.edu/static/models/iML1515.json", "iML1515.json");

using COBREXA
using SparseArrays
using LinearAlgebra
using Tulip

model = load_model(CoreModel, "iML1515.json")

# ## Describe the variants
# Each variant is described by the substrate it consumes for energy/carbon

## Each variant is specified by the substrate it consumes for energy/carbon
variants = [
    "glc__D_e",
    "xyl__D_e",
    "fru_e",
    "man_e",
    "gal_e",
    "tre_e",
    "malt_e",
    "lcts_e",
    "fuc__L_e",
    "arab__L_e",
    "ac_e",
    "etoh_e",
    "lac__D_e",
]

n_models = length(variants) # number of different models to construct and merge

# ## Build the community stoichiometric matrix
# The community stoichiometric matrix is constructed by combining the
# stoichiometric matrices of the individual models, as well as adding
# environmental variables.

ex_rxn_inds = find_exchange_reactions(model; exclude_biomass = true)
ex_met_inds = find_exchange_metabolites(model; exclude_biomass = true)
n_env_vars = length(ex_rxn_inds) # number of exchange metabolites/reactions are equal
n_rows, n_cols = size(stoichiometry(model))

S = spzeros(n_models * n_rows + n_env_vars, n_models * n_cols + n_env_vars)

## put individual model stoichiometric matrices along diagonals
for n = 1:n_models
    row_start = (n - 1) * n_rows + 1
    row_end = n_rows + (n - 1) * n_rows
    col_start = (n - 1) * (n_cols) + 1
    col_end = n_cols + (n - 1) * (n_cols)
    S[row_start:row_end, col_start:col_end] .= stoichiometry(model)
end

## add environmental exchanges
S[n_models*n_rows+1:end, n_models*n_cols+1:end] .= -sparse(I, n_env_vars, n_env_vars) # ENV met -> ∅

## modify individual exchange reactions
for n = 1:n_models
    rows = collect(1:n_env_vars) .+ n_models * n_rows
    cols = (n - 1) * n_cols .+ ex_rxn_inds
    for (i, j) in zip(rows, cols)
        S[i, j] = 1.0 # then EXT met -> ENV met
    end
end

# ## Build flux bound vectors
# Each flux needs to be constrained. For the most part, use the flux constraints
# contained in the default model.
lbs = spzeros(n_models * n_cols + n_env_vars)
ubs = spzeros(n_models * n_cols + n_env_vars)
single_lbs, single_ubs = bounds(model) # same bounds for each organism
for i = 1:n_models
    col_start = (i - 1) * (n_cols) + 1
    col_end = n_cols + (i - 1) * (n_cols)
    lbs[col_start:col_end] .= single_lbs
    ubs[col_start:col_end] .= single_ubs
end

## match exchange bounds to environmental fluxes for simplicity
lbs[n_models*n_cols+1:end] .= single_lbs[ex_rxn_inds]
ubs[n_models*n_cols+1:end] .= single_ubs[ex_rxn_inds];

# ## Create new names for the model's reactions and metabolites
# Append information about the variant to the reactions and metabolites of each
# model. This makes it easier to inspect/understand the flux solution at the
# end.

env_mets = metabolites(model)[ex_met_inds]
## name reactions
rxn_ids = [
    vcat([reactions(model) .* "_org_$(variants[i])" for i = 1:n_models]...)
    "EX_" .* env_mets .* "_ENV"
]
#
## name metabolites
met_ids = [
    vcat([metabolites(model) .* "_org_$(variants[i])" for i = 1:n_models]...)
    env_mets .* "_ENV"
]

# ##  Adjust variants 
# Each variant can only consume the substrate associated with it for
# carbon/energy. Adjust the constraints to reflect this.

for n = 1:n_models

    ## can import only the variant specific substrate
    ind = first(indexin(["EX_$(variants[n])_org_$(variants[n])"], rxn_ids))
    lbs[ind] = -1000.0
    ubs[ind] = 0.0

    ## set environmental availability of sugar (by default these fluxes are export only)
    if !(variants[n] in ["ac_e", "etoh_e", "lac__D_e"])
        ind = first(indexin(["EX_$(variants[n])_ENV"], rxn_ids))
        lbs[ind] = -10.0
        ubs[ind] = 0.0
    else # set import bound for non-sugar substrates
        ind = first(indexin(["EX_$(variants[n])_org_$(variants[n])"], rxn_ids))
        lbs[ind] = -50.0
        ubs[ind] = 0.0
    end

    if n != 1 # remove ability to metabolize glucose, the default condition of the model
        ind = first(indexin(["EX_$(variants[1])_org_$(variants[n])"], rxn_ids))
        lbs[ind] = 0.0
        ubs[ind] = 0.0
    end
end

# ## Add biomass objective function
# A biomass reaction is added that ensures that each microbe will grow at the
# same rate. In short, this is achieved by creating new metabolites, which are
# just the biomass of each model. These biomass metabolites are then combined to
# create an overall biomass reaction. Since the proportions are locked into the
# reaction (weighted), each microbe must grow at the same rate.

## create biomass metabolites
S = vcat(S, [spzeros(size(S, 2))' for n = 1:n_models]...)
met_ids = [met_ids; ["biomass_org_$(variants[n])" for n = 1:n_models]]
for n = 1:n_models
    row_num = first(indexin(["biomass_org_$(variants[n])"], met_ids))
    col_num = first(indexin(["BIOMASS_Ec_iML1515_core_75p37M_org_$(variants[n])"], rxn_ids))
    S[row_num, col_num] = 1.0 # creates 1 biomass
end
## create community biomass objective function
obj_rxn = spzeros(size(S, 1))
biomass_mets = findall(x -> occursin("biomass_org", x), met_ids) # a WT bof also exists
obj_rxn[biomass_mets] .= -1.0
S = [S obj_rxn]
rxn_ids = [rxn_ids; "community_biomass"]

# ## Create community model
# A community model is created by combining the elements of the model created
# thus far. `CoreModel` is used since it is made for large-scale community
# analysis.

lbs = [lbs; 0.0]
ubs = [ubs; 1000.0]
c = spzeros(length(lbs))
c[end] = 1.0
b = spzeros(size(S, 1))
dropzeros!(lbs)
dropzeros!(ubs)
dropzeros!(S)

community_model = CoreModel(S, b, c, lbs, ubs, rxn_ids, met_ids)

# ## Analyze community model
d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)

bof_rxn_inds = findall(x -> occursin("BIOMASS_Ec_iML1515_core", x), rxn_ids)
for obj_id in rxn_ids[bof_rxn_inds]
    println(obj_id, ": ", d[obj_id])
end

# 
for n = 1:n_models
    env_ex = "EX_" * variants[n] * "_org_$(variants[n])"
    println(env_ex, ": ", d[env_ex])
end
