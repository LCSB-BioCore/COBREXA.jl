# # Community FBA models

using COBREXA

# Here we will construct a community FBA model of two  *E. coli* "core" models
# that can interact by exchanging selected metabolites. To do this, we will need
# the model, which we can download if it is not already present.

import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use Tulip here:

import JSONFBCModels
import Tulip
import AbstractFBCModels as A
import ConstraintTrees as C

model = load_model("e_coli_core.json")

# Community models work by joining its members together through their exchange
# reactions, weighted by the abundance of each microbe. These exchange reactions
# are then linked to an environmental exchange. For more theoretical details,
# see "Gottstein, et al, 2016, Constraint-based stoichiometric modelling from
# single organisms to microbial communities, Journal of the Royal Society
# Interface". 

# ## Building a community of two *E. coli*s 

# Here we will construct a simple community of two interacting microbes. To do
# this, we need to import the models. We import the models are ConstraintTrees,
# because it is easier to build the model explicitly than rely on an opaque
# one-shot function.

ecoli1 = fbc_model_constraints(model)
ecoli2 = fbc_model_constraints(model)

# Since the models are joined through their individual exchange reactions to an
# environmental exchange reactionq, we need to identify all possible exchange
# reactions in the community. Since the models are the same, this is
# straightforward here. Additionally, we need to specify the upper and lower
# bounds of these environmental exchange reactions.
lbs, ubs = A.bounds(model)

env_ex_rxns = Dict(
    rid => (lbs[i], ubs[i]) for
    (i, rid) in enumerate(A.reactions(model)) if startswith(rid, "EX_")
)

# Now we simply create an blank model that only includes environmental exchange reactions.

m = build_community_environment(env_ex_rxns)

# Next we join each member microbe to the model.
m += :bug1^ecoli1
m += :bug2^ecoli2

# We also need to specify the abundances of each member, as this weights the
# flux of each metabolite each member microbe can share with other members or
# the environment.
member_abundances = [(:bug1, 0.2), (:bug2, 0.8)]

m *= :environmental_exchange_balances^link_environmental_exchanges(m, member_abundances)

# Finally, the most sensible community FBA simulation involves assuming the
# growth rate of the models is the same. In this case, we simply set the growth
# rate flux of each member to be the same.
m *=
    :equal_growth_rate_constraint^equal_growth_rate_constraints([
        (:bug1, m.bug1.fluxes.:BIOMASS_Ecoli_core_w_GAM.value),
        (:bug2, m.bug2.fluxes.:BIOMASS_Ecoli_core_w_GAM.value),
    ])

# Since each growth rate is the same, we can pick any of the growth rates as the
# objective for the simulation.
m *= :objective^C.Constraint(m.bug1.fluxes.:BIOMASS_Ecoli_core_w_GAM.value)

# Since the models are usually used in a mono-culture context, the glucose input
# for each individual member is limited. We need to undo this limitation, and
# rather rely on the constrained environmental exchange reaction (and the bounds
# we set for it earlier).
m.bug1.fluxes.EX_glc__D_e.bound = (-1000.0, 1000.0)
m.bug2.fluxes.EX_glc__D_e.bound = (-1000.0, 1000.0)

# We can also be interesting, and limit respiration in one of the members, to
# see what effect this has on the community.
m.bug1.fluxes.CYTBD.bound = (-10.0, 10.0)

# Finally, we can simulate the system!
sol = optimized_constraints(
    m;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [set_optimizer_attribute("IPM_IterationsLimit", 1000)],
)

@test isapprox(sol.:objective, 0.66686196344, atol = TEST_TOLERANCE) #src

# At the moment the members cannot really exchange any metabolites. We can
# change this by changing their individual exchange bounds. 
mets = [:EX_akg_e, :EX_succ_e, :EX_pyr_e, :EX_acald_e, :EX_fum_e, :EX_mal__L_e]
for met in mets
    m.bug1.fluxes[met].bound = (-1000.0, 1000.0)
    m.bug2.fluxes[met].bound = (-1000.0, 1000.0)
end

sol = optimized_constraints(
    m;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [set_optimizer_attribute("IPM_IterationsLimit", 1000)],
)


# We can see that by allowing the microbes to share metabolites, the growth rate
# of the system as a whole increased! We can inspect the individual exchanges to
# see which metabolites are being shared (pyruvate in this case).
bug1_ex_fluxes = Dict(k => v for (k, v) in sol.bug1.fluxes if startswith(string(k), "EX_"))
bug2_ex_fluxes = Dict(k => v for (k, v) in sol.bug2.fluxes if startswith(string(k), "EX_"))

#!!! warning "Flux units"
# The unit of the environmental exchange reactions (mmol/gDW_total_biomass/h) is
# different to the unit of the individual species fluxes
# (mmol/gDW_species_biomass/h). This is because the mass balance needs to take
# into account the abundance of each species for the simulation to make sense.
# In this specific case, look at the flux of pyruvate (EX_pyr_e). There is no
# environmental exchange flux, so the two microbes share the metabolite.
# However,  `bug1_ex_fluxes[:EX_pyr_e] != bug2_ex_fluxes[:EX_pyr_e]`, but rather
# `abundance_bug1 * bug1_ex_fluxes[:EX_pyr_e] != abundance_bug2 *
# bug2_ex_fluxes[:EX_pyr_e]`. Take care of this when comparing fluxes!

@test isapprox(
    abs(0.2 * bug1_ex_fluxes[:EX_pyr_e] + 0.8 * bug2_ex_fluxes[:EX_pyr_e]),
    0.0,
    atol = TEST_TOLERANCE,
) #src
