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

variants = ["EX_glc__D_e", "EX_xyl__D_e"]
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

models = [load_model(CoreModel, "iML1515.json") for i = 1:n_models];

ex_rxn_inds = find_exchange_reactions(models[1]; exclude_biomass = true)
ex_met_inds = find_exchange_metabolites(models[1]; exclude_biomass = true)
n_env_vars = length(ex_rxn_inds) # number of exchange metabolites/reactions are equal
n_rows, n_cols = size(stoichiometry(models[1]))

# create the community stoichiometric matrix
community_S = spzeros(n_models * n_rows + n_env_vars, n_models * n_cols + n_env_vars)

# place along diagonals
for n = 1:n_models
    row_start = (n - 1) * n_rows + 1
    row_end = n_rows + (n - 1) * n_rows
    col_start = (n - 1) * (n_cols) + 1
    col_end = n_cols + (n - 1) * (n_cols)
    community_S[row_start:row_end, col_start:col_end] .= stoichiometry(models[n])
end

# add environmental exchanges
community_S[n_models*n_rows+1:end, n_models*n_cols+1:end] .=
    -sparse(I, n_env_vars, n_env_vars) # ENV met -> ∅

# modify individual exchange reactions
for n = 1:n_models
    rows = collect(1:n_env_vars) .+ n_models * n_rows
    cols = (n - 1) * n_cols .+ ex_rxn_inds
    for (i, j) in zip(rows, cols)
        community_S[i, j] = 1.0 # then EXT met -> ENV met
    end
end
community_S

# build flux bound vectors
lbs = spzeros(n_models*n_cols+n_env_vars)
ubs = spzeros(n_models*n_cols+n_env_vars)
single_lbs, single_ubs = bounds(models[1]) # same bounds for each organism
for i = 1:n_models
    col_start = (i - 1) * (n_cols) + 1
    col_end = n_cols + (i - 1) * (n_cols)
    lbs[col_start:col_end] .= single_lbs
    ubs[col_start:col_end] .= single_ubs
end

# match exchange bounds
lbs[n_models*n_cols+1:end] .= single_lbs[ex_rxn_inds]
ubs[n_models*n_cols+1:end] .= single_ubs[ex_rxn_inds]

env_mets = metabolites(models[1])[ex_met_inds]
# name reactions
rxn_ids = [vcat([reactions(models[i]).*"_org$i" for i=1:n_models]...); "EX_".*env_mets.*"_ENV"]
# name metabolites
met_ids = [vcat([reactions(models[i]).*"_org$i" for i=1:n_models]...); env_mets.*"_ENV"]

# adjust for variants
for n in 1:n_models

    ind = first(indexin(["$(variants[n])_org$n"], rxn_ids))
    lbs[ind] = -1000.0 # can only import this metabolite
    ubs[ind] = 0.0    

    ind = first(indexin(["$(variants[n])_ENV"], rxn_ids))
    lbs[ind] = -10.0
    ubs[ind] = 0.0    
    
    if n != 1 # remove ability to metabolize glucose
        ind = first(indexin(["EX_glc__D_e_org$n"], rxn_ids))
        lbs[ind] = 0.0
        ubs[ind] = 0.0
    end
end

# biomass objective function
c = spzeros(length(lbs))
obj_func_inds
c[obj_func_inds] = 1.0 

dropzeros!(lbs)
dropzeros!(ubs)
dropzeros!(community_S)
