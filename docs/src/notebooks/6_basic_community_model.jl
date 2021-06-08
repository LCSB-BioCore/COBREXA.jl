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

# ## Create mutants by knocking out certain reactions

base_model = load_model(CoreModel, "e_coli_core.json")

cytbd_knockout = remove_reactions(base_model, "CYTBD")
sol = flux_balance_analysis_dict(cytbd_knockout, Tulip.Optimizer)
sol["BIOMASS_Ecoli_core_w_GAM"]

atps4r_knockout = remove_reactions(base_model, "ATPS4r")
sol = flux_balance_analysis_dict(atps4r_knockout, Tulip.Optimizer)
sol["BIOMASS_Ecoli_core_w_GAM"]

pdh_knockout = remove_reactions(base_model, "PDH")
sol = flux_balance_analysis_dict(pdh_knockout, Tulip.Optimizer)
sol["BIOMASS_Ecoli_core_w_GAM"]

ex_rxns = filter(looks_like_exchange_reaction, reactions(base_model))
ex_mets = [first(keys(reaction_stoichiometry(base_model, ex_rxn))) for ex_rxn in ex_rxns]

models = [cytbd_knockout, atps4r_knockout, pdh_knockout]
community_model = join_with_exchanges(models, ex_rxns, ex_mets; 
        add_biomass_objective=true, 
        biomass_ids=["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"], 
        model_names=["cytbd_ko", "atps4r_ko", ]
    ) 
