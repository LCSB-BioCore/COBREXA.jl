# # Growth media analysis

# Nutrient availability is a major driving factor for growth of microorganisms
# and energy production in cells. Here, we demonstrate two main ways to examine
# the nutrient consumption with COBREXA.jl: Simulating deficiency of nutrients,
# and finding the minimal flux of nutrients required to support certain model
# output.

# As always, we work on the toy model of *E. coli*:

using COBREXA, GLPK

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

model = load_model(ObjectModel, "e_coli_core.xml")

# ## What nutrients does my model need to grow?

# The models usually ingest nutrients through exchange reactions. By changing
# the bounds on the exchange reactions, you can limit the intake of the
# nutrients and thus simulate the nutrient deficiency. If applied
# programatically to multiple exchanges, this can give you a good overview of
# what nutrients impact the model most.
#
# To check the viability of a single nutrient, you can simply change a bound on
# a selected exchange reaction and simulate the model with a limited amount.

biomass = "R_BIOMASS_Ecoli_core_w_GAM"

model_limited = change_bound(model, "R_EX_glc__D_e", lower_bound = -1.0)

#md # !!! tip "Exchange directions"
#md #     By a convention, the direction of exchange reaction usually goes from the
#md #     model into the environment, representing the "production". Limiting the
#md #     intake thus happens by disabling the "negative production", i.e., placing
#md #     a lower bound.

original_production = flux_balance_analysis_dict(model, GLPK.Optimizer)[biomass]
limited_production = flux_balance_analysis_dict(model_limited, GLPK.Optimizer)[biomass]

original_production, limited_production

# Function [`flux_summary`](@ref) can help with quickly spotting what has
# changed:

flux_summary(flux_balance_analysis_dict(model_limited, GLPK.Optimizer))

# Similarly, you can check that the model can survive without oxygen, at the cost
# of switching the metabolism to ethanol production:

flux_summary(
    flux_balance_analysis_dict(
        change_bound(model, "R_EX_o2_e", lower_bound = 0.0),
        GLPK.Optimizer,
    ),
)

# The effect of all nutrients on the metabolism can be scanned using [`screen`](@ref). The [`change_bound`](@ref) function is, for this purpose, packed in a variant specified [`with_changed_bound`](@ref):

exchanges = filter(looks_like_exchange_reaction, variables(model))

exchanges .=> screen(
    model,
    variants = [
        [with_changed_bound(exchange, lower_bound = 0.0)] for exchange in exchanges
    ],
    analysis = m -> begin
        res = flux_balance_analysis_dict(m, GLPK.Optimizer)
        isnothing(res) ? nothing : res[biomass]
    end,
)


# Similarly to gene knockouts, you can also examine the effect of combined
# nutrient deficiencies. To obtain a more interesting result, we may examine
# the effect of slight deficiencies of pairs of intake metabolites. For
# simplicity, we show the result only on a small subset of the exchanges:

selected_exchanges = [
    "R_EX_pi_e",
    "R_EX_gln__L_e",
    "R_EX_nh4_e",
    "R_EX_pyr_e",
    "R_EX_fru_e",
    "R_EX_glu__L_e",
    "R_EX_glc__D_e",
    "R_EX_o2_e",
]

screen(
    model,
    variants = [
        [with_changed_bounds([e1, e2], lower_bound = [-1.0, -0.1])] for
        e1 in selected_exchanges, e2 in selected_exchanges
    ],
    analysis = m -> begin
        res = flux_balance_analysis_dict(m, GLPK.Optimizer)
        isnothing(res) ? nothing : res[biomass]
    end,
)

# The result describes combinations of nutrient deficiencies -- the nutrient
# that corresponds to the row is mildly deficient (limited to uptake 1.0), and
# the one that corresponds to the column is severely limited (to uptake 0.1).

#md # !!! tip "Screening can be easily parallelized"
#md #     To speed up larger analyses, remember that execution of [`screen`](@ref)
#md #     can be [parallelized to gain speedup](../distributed/1_functions.md). Parallelization in `screen` is optimized to avoid
#md #     unnecessary data transfers that may occur when using trivial `pmap`.

# ## What is the minimal flux of nutrients for my model to grow?

# You can compute the minimal flux (i.e., mass per time) of required nutrients
# by constraining the model growth to a desired lower bound, and then optimize
# the model with an objective that minimizes intake of all exchanges (i.e.,
# given the directionality convention of the exchanges, actually maximizes the
# flux through all exchange reactions along their direction).

model_with_bounded_production = change_bound(model, biomass, lower_bound = 0.1) #minimum required growth

minimal_intake_production = flux_balance_analysis_dict(
    model_with_bounded_production,
    GLPK.Optimizer,
    modifications = [change_objective(exchanges)],
);

# Metabolite "cost" data may be supplemented using the `weights` argument of
# [`change_objective`](@ref), to reflect e.g. the different molar masses or
# energetic values of different nutrients.
#
# In our simple case, we obtain the following minimal required intake:

flux_summary(minimal_intake_production)
