# # Building and analysing a small community model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we will use `COBREXA` to build and analyze a small community model
# consisting of three *E. coli* mutants using the `CoreModel`. We will use an
# objective function that enforces equal growth rates.

# We will first construct a community of only two mutants to illustrate the
# effect of the community biomass objective function. Then we will add a third
# member that has a lethal knockout. Due to the bounds on the exchange reactions
# these three models are incapable of sharing resources - hence the maximum
# growth rate will be zero. By changing the bounds we can allow resource sharing,
# "saving" the community.

# ## Load the base model

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

using COBREXA
using Tulip

# ## Load the models and inspect fba solutions

base_model = load_model(CoreModel, "e_coli_core.json") # base from from which the knockouts will be constructed

cytbd_knockout_model = remove_reaction(base_model, "CYTBD") # knockout the CYTBD (cytochrome oxidase) reaction
sol = flux_balance_analysis_dict(cytbd_knockout_model, Tulip.Optimizer)
sol["BIOMASS_Ecoli_core_w_GAM"] # Cytochrome oxidase knockout μ (growth rate)
#
atps4r_knockout_model = remove_reaction(base_model, "ATPS4r") # knockout the ATP synthase reaction
sol = flux_balance_analysis_dict(atps4r_knockout_model, Tulip.Optimizer)
sol["BIOMASS_Ecoli_core_w_GAM"] # ATP synthase knockout μ
#
eno_knockout_model = remove_reaction(base_model, "ENO") # knockout the enolase reaction
sol = flux_balance_analysis_dict(eno_knockout_model, Tulip.Optimizer)
sol["BIOMASS_Ecoli_core_w_GAM"] # Enolase knockout μ, cannot grow by itself

# ## Build a community model of the cytochrome oxidase knockout and the ATP synthase knockout models

ex_rxn_mets = Dict(
    ex_rxn => first(keys(reaction_stoichiometry(base_model, ex_rxn))) for
    ex_rxn in filter(looks_like_exchange_reaction, reactions(base_model))
) # identify exchange reactions heuristically
#
model_names = ["cytbd_ko", "atps4r_ko"]
community_model = join_with_exchanges(
    CoreModel,
    [cytbd_knockout_model, atps4r_knockout_model],
    ex_rxn_mets;
    biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"],
    model_names = model_names,
)

# ## Set exchange reaction bounds of community model based on the bounds of the individual models

env_ex_rxn_idxs = indexin(keys(ex_rxn_mets), reactions(community_model)) # identify the global (environmental exchange reactions)
cytbd_ex_rxn_idxs = indexin(keys(ex_rxn_mets), reactions(cytbd_knockout_model)) # identify the indices of the corresponding exchange reactions in the original models
atps4r_ex_rxn_idxs = indexin(keys(ex_rxn_mets), reactions(atps4r_knockout_model))

# In case some exchange reactions are not present in both models, set
# environmental exchange bound to the sum of the individual exchange bounds
for (env_ex, m2_ex, m1_ex) in zip(env_ex_rxn_idxs, cytbd_ex_rxn_idxs, atps4r_ex_rxn_idxs)
    m2lb = isnothing(m2_ex) ? 0.0 : atps4r_knockout_model.xl[m2_ex]
    m2ub = isnothing(m2_ex) ? 0.0 : atps4r_knockout_model.xu[m2_ex]
    m1lb = isnothing(m1_ex) ? 0.0 : cytbd_knockout_model.xl[m1_ex]
    m1ub = isnothing(m1_ex) ? 0.0 : cytbd_knockout_model.xu[m1_ex]
    change_bounds!(community_model, [env_ex]; lower = [m1lb + m2lb], upper = [m1ub + m2ub])
end

# ## Add objective function to community model`

biomass_ids = model_names .* "_BIOMASS_Ecoli_core_w_GAM"
update_objective!(community_model, biomass_ids, objective_id = "community_biomass")

# ## Perform community FBA

d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)
d["community_biomass"] # community μ
# Notice, the growth rate is limited to the slowest organism as per the objective function

# ## Add the enolase knockout to the community model

community_model = add_model_with_exchanges(
    community_model,
    eno_knockout_model,
    ex_rxn_mets;
    model_name = "eno_ko",
    biomass_id = "BIOMASS_Ecoli_core_w_GAM",
)

push!(model_names, "eno_ko")
biomass_ids = model_names .* "_BIOMASS_Ecoli_core_w_GAM"
update_objective!(community_model, biomass_ids, "community_biomass")

d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)
println("Community μ = ", d["community_biomass"])
# Notice that the high communal growth rate is 0, due to the enolase knockout.
# The reason for this behaviour: enolase is a central reaction in glycolysis - without
# it the organism cannot access the lower glycolysis pathways or the TCA cycle, hence
# the model predicts no growth for the knockout, and hence no growth for the system since
# they all have to have the same growth rate.

# ## Allow the mutants to rescue each other by sharing pyruvate

pyr_exs = model_names .* "_EX_pyr_e"
change_bounds!(community_model, pyr_exs; lower = fill(-1000.0, 3), upper = fill(1000.0, 3))

d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)
d["community_biomass"] # community μ
# Notice that the growth rate is now above 0! Nutrient sharing saved the day!
