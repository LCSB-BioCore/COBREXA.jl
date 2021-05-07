# # Building and analysing a community model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this tutorial we will use `COBREXA` to build and analyze a community model  
# consisting of multiple variants of *E. coli* knockouts using the `CoreModel`.
# Here each knockout will only be able to metabolize one sugar. 

# If it is not already present, load a large scale *E. coli* model.

!isfile("iML1515.xml") &&
    download("http://bigg.ucsd.edu/static/models/iML1515.xml", "iML1515.xml")
#
using COBREXA
using SparseArrays
using LinearAlgebra
using Tulip 

variants = ["EX_glc__D_e",
            "EX_xyl__D_e",]
            # "EX_fru_e", 
            # "EX_sucr_e",
            # "EX_man_e",
            # "EX_gal_e",
            # "EX_tre_e",
            # "EX_malt_e",
            # "EX_lcts_e",
            # "EX_fuc__L_e",
            # "EX_arab__L_e"] # each variant can only metabolize the sugar it is associated with

n_models = length(variants)

models = [load_model(CoreModel, "iML1515.mat") for i=1:n_models];

# since each model is identical, all the exchange reactions will be 
# in the same indices
ex_rxn_inds = find_exchange_reactions(models[1]; exclude_biomass = true) 
ex_met_inds = find_exchange_metabolites(models[1]; exclude_biomass = true)
n_env_vars = length(ex_rxn_inds) # number of exchange metabolites/reactions are equal
n_rows, n_cols = size(stoichiometry(models[1]))

# create the community stoichiometric matrix
community_S = spzeros(n_models*n_rows + n_env_vars, n_models*n_cols + n_env_vars)

# place along diagonals
for i=1:n_models
    row_start = (i-1)*n_rows + 1
    row_end = n_rows + (i-1)*n_rows
    col_start = (i-1)*(n_cols) + 1
    col_end =  n_cols + (i-1)*(n_cols)
    community_S[row_start:row_end, col_start:col_end] .= stoichiometry(models[i])
end

# add environmental exchanges
community_S[n_models*n_rows+1:end,n_models*n_cols+1:end] .= -sparse(I, n_env_vars, n_env_vars) # ENV met -> ∅

# modify individual exchange reactions
for n=1:n_models
    rows = collect(1:n_env_vars) .+ n_models*n_rows
    cols = (n-1)*n_cols .+ ex_rxn_inds
    for (i, j) in zip(rows, cols)
        community_S[i, j] = 1.0 # then EXT met -> ENV met
    end
end
community_S

single_lbs, single_ubs = bounds(models[1])
lbs = sparse([])
ubs = spzeros(n_cols(n_models))

# constrain each individual model
for i=1:n_models
    # allow only variant sugar as carbon source
    ind = first(indexin([variants[i]], reactions(models[i])))
    models[i].xl[ind] = -1000.0
    models[i].xu[ind] = 0.0

    if i != 1 # ensure glucose can't be used as a carbon source for the other bugs
        ind = first(indexin([variants[1]], reactions(models[i])))
        models[i].xl[ind] = 0.0
        models[i].xu[ind] = 0.0    
    end
end

# create environmental metabolites

first(indexin([variants[i]], reactions(models[i])))
i = 1
d = flux_balance_analysis_dict(models[i], Tulip.Optimizer)
d["R_BIOMASS_Ec_iML1515_core_75p37M"]
d[variants[i]]
d[variants[1]]

# m = load_model("iML1515.xml") 
m = load_model(joinpath("..","models","e_coli_core.xml")) 
d = flux_balance_analysis_dict(m, Tulip.Optimizer)
bounds(m)
