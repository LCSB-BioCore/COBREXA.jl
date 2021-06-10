# # Building and analysing a small community model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we will use `COBREXA` to build and analyze a small community model  
# consisting of three *E. coli* mutants using the `CoreModel`. We will use an
# objective function that enforces equal growth rates.

# ## Load the base model

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

using COBREXA
using Tulip

# ## Load the models and inspect fba solutions

base_model = load_model(CoreModel, "e_coli_core.json")

cytbd_knockout = remove_reactions(base_model, "CYTBD")
sol = flux_balance_analysis_dict(cytbd_knockout, Tulip.Optimizer)
println("Cytochrome oxidase knockout μ = ", sol["BIOMASS_Ecoli_core_w_GAM"])

atps4r_knockout = remove_reactions(base_model, "ATPS4r")
sol = flux_balance_analysis_dict(atps4r_knockout, Tulip.Optimizer)
println("ATP synthase knockout μ = ", sol["BIOMASS_Ecoli_core_w_GAM"])

eno_knockout = remove_reactions(base_model, "ENO")
sol = flux_balance_analysis_dict(eno_knockout, Tulip.Optimizer)
println("Enolase knockout μ = ", sol["BIOMASS_Ecoli_core_w_GAM"])
# notice that the enolase mutant cannot grow by itself.

# ## Build a community model of the Cytochrome oxidase and ATP synthase knockouts 

ex_rxns = filter(looks_like_exchange_reaction, reactions(base_model))
ex_mets = [first(keys(reaction_stoichiometry(base_model, ex_rxn))) for ex_rxn in ex_rxns]

models = [cytbd_knockout, atps4r_knockout]
community_model = join_with_exchanges(models, ex_rxns, ex_mets; 
        add_biomass_objective=true, 
        biomass_ids=["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"], 
        model_names=["cytbd_ko", "atps4r_ko"]
        )     

# ## Set bounds of community model

env_ex_rxn_idxs = indexin(ex_rxns, reactions(community_model))
cytbd_ex_rxn_idxs = indexin(ex_rxns, reactions(cytbd_knockout))
atps4r_ex_rxn_idxs = indexin(ex_rxns, reactions(atps4r_knockout))

for (env_ex, m2_ex, m1_ex) in zip(env_ex_rxn_idxs, cytbd_ex_rxn_idxs, atps4r_ex_rxn_idxs)
    # in case some exchange reactions are not present in both models
    m2lb = isnothing(m2_ex) ? 0.0 : atps4r_knockout.xl[m2_ex]
    m2ub = isnothing(m2_ex) ? 0.0 : atps4r_knockout.xu[m2_ex]

    m1lb = isnothing(m1_ex) ? 0.0 : cytbd_knockout.xl[m1_ex]
    m1ub = isnothing(m1_ex) ? 0.0 : cytbd_knockout.xu[m1_ex]

    change_bounds!(community_model, [env_ex]; xl= [m1lb + m2lb], xu = [m1ub + m2ub])
end

# ## Add objective function to community model

biomass_ids = ["cytbd_ko_BIOMASS_Ecoli_core_w_GAM", "atps4r_ko_BIOMASS_Ecoli_core_w_GAM"]
add_objective!(
    community_model,
    biomass_ids;
    objective_column_index = first(
        indexin(["community_biomass"], reactions(community_model)),
    ),
)

# ## Perform community FBA

d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)
println("Community μ = ", d["community_biomass"])
# The growth rate is limited to the slowest organism as per the objective function

# ## Add the Enolase knockout to the community model

community_model = add_model_with_exchanges(
    community_model,
    eno_knockout,
    ex_rxns,
    ex_mets;
    model_name = "eno_ko",
    biomass_id = "BIOMASS_Ecoli_core_w_GAM",
)

# ## update the biomass reaction to include the new mutant

biomass_ids = ["cytbd_ko_BIOMASS_Ecoli_core_w_GAM", "atps4r_ko_BIOMASS_Ecoli_core_w_GAM", "eno_ko_BIOMASS_Ecoli_core_w_GAM"]
add_objective!(
    community_model,
    biomass_ids;
    objective_column_index = first(
        indexin(["community_biomass"], reactions(community_model)),
    ),
)

d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)
println("Community μ = ", d["community_biomass"])   
# Notice that the high communal growth rate is 0, due to the enolase knockout.

# ## Allow the mutants to share pyruvate and rescue each other

pyr_exs = ["cytbd_ko_EX_pyr_e", "atps4r_ko_EX_pyr_e", "eno_ko_EX_pyr_e"]
change_bounds!(community_model, pyr_exs; xl = repeat([-1000.0], inner=3), xu = repeat([1000.0], inner=3))

d = flux_balance_analysis_dict(
    community_model,
    Tulip.Optimizer;
    modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
)
println("Community μ = ", d["community_biomass"])   
