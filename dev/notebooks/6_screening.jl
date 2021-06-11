# # Exploring model variants with `screen`

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# This notebooks demonstrates a simple use of [`screen`](@ref) to explore large
# number of model variants. On the toy *E. Coli* model, we try to map the
# impact of knocking out single reactions and 2-reaction combinations.

# First, let's download the data and load the packages and the model
!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

using COBREXA, Tulip

model = load_model(StandardModel, "e_coli_core.json")

# ## Preparing the functions
#
# While we could use the [`with_set_bound`](@ref) to limit the reaction rates,
# but we will make a slightly more precise, usage-tailored modification. This
# is a straightforward modification of the [`with_set_bound`](@ref) that does
# not set bounds "outside" of the original bounds:

with_limited_rate(reaction_id::String, limit) =
    model::StandardModel -> begin
        m = copy(model)
        m.reactions = copy(model.reactions)
        r = m.reactions[reaction_id] = copy(model.reactions[reaction_id])
        if -limit > r.lb
            r.lb = -limit
        end
        if limit < r.ub
            r.ub = limit
        end
        m
    end

# ## Knocking out single reactions
#
# This can be applied to see how much biomass can the model produce with
# certain reactions limited to almost zero:

get_biomass(x) = isnothing(x) ? 0 : x["BIOMASS_Ecoli_core_w_GAM"]

productions = screen_variants(
    model,
    [[with_limited_rate(rxn, 0.1)] for rxn in reactions(model)],
    model -> get_biomass(flux_balance_analysis_dict(model, Tulip.Optimizer)),
)

# This can be nicely plotted to give a more comprehensive overview of which
# reactions are critical or not:

using Plots
bar(reactions(model), productions, orientation = :hor, dpi = 600)

# ## Knocking out reaction combinations
#
# It is very easy to prepare matrices of biomass productions from all possible
# two-reaction knockouts. To make it more interesting, we will restrict one of
# the reactions of the pair a bit less, to see more possible outcomes.

# We do not process all reactions here to make the notebook rendering
# efficient, but you can easily remove the restriction, and speed the process
# up using parallel execution, by specifying `workers` argument (see
# documentation of [`screen`](@ref) for details)

rxns = reactions(model)

productions = screen_variants(
    model,
    [
        [with_limited_rate(rxn1, 3), with_limited_rate(rxn2, 0.1)] for rxn1 in rxns,
        rxn2 in rxns
    ],
    model -> get_biomass(flux_balance_analysis_dict(model, Tulip.Optimizer)),
)

heatmap(productions, dpi = 600)
