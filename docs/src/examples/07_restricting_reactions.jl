
# # Restricting and disabling individual reactions

# Here, we show several methods how to explore the effect of disabling or choking the reactions in the models.

# First, download the demonstration data and load the packages as usual:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model(ObjectModel, "e_coli_core.json")

# ## Disabling a reaction

# There are several possible ways to disable a certain reaction in the model.
# The easiest way is to use [`change_bound`](@ref) or [`change_bounds`](@ref)
# to create a variant of the model that has the corresponding bounds modified
# (or, alternatively, a pipeable "variant" version
# [`with_changed_bound`](@ref)).
#
# Alternatively, you could utilize [`modify_constraint`](@ref) as a
# modification that acts directly on the JuMP optimization model. That may be
# useful if you first apply some kind of complicated constraint scheme
# modification, such as [`add_loopless_constraints`](@ref).
#
# In turn, in the simple case, the following 2 ways of disabling the FBA
# reaction are equivalent. The first, making a variant of the model structure,
# might be slightly preferred because it better composes with other changes;
# the second does not compose as well but may be more efficient (and thus
# faster) in certain situations:

flux1 = flux_balance_analysis_vec(
    model |> with_changed_bound("FBA", lower_bound = 0.0, upper_bound = 0.0),
    GLPK.Optimizer,
);

flux2 = flux_balance_analysis_vec(
    model,
    GLPK.Optimizer,
    modifications = [modify_constraint("FBA", lb = 0.0, ub = 0.0)],
);

# The solutions should not differ a lot:
sum((flux1 .- flux2) .^ 2)

# ## Restricting a reaction
#
# Quite naturally, you can restruct the reaction to a limited flow, simulating
# e.g. nutrient deficiency:

original_flux = flux_balance_analysis_dict(model, GLPK.Optimizer);

restricted_flux = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer,
    modifications = [modify_constraint("EX_o2_e", lb = -0.1, ub = 0.0)],
);

# The growth in the restricted case is, expectably, lower than the original one:
original_flux["BIOMASS_Ecoli_core_w_GAM"], restricted_flux["BIOMASS_Ecoli_core_w_GAM"]

# ## Screening for sensitive and critical reactions
#
# Using higher-order analysis scheduling functions (in particular
# [`screen`](@ref)), you can easily determine which reactions play a crucial
# role for the model viability and which are not very important.

# We can take all reactions where the flux is not zero:

running_reactions = [(rid, x) for (rid, x) in original_flux if abs(x) > 1e-3]

# ...and choke these reactions to half that flux, computing the relative loss
# of the biomass production::

screen(
    model,
    variants = [
        [with_changed_bound(rid, lower_bound = -0.5 * abs(x), upper_bound = 0.5 * abs(x))]
        for (rid, x) in running_reactions
    ],
    args = running_reactions,
    analysis = (m, rid, _) ->
        rid =>
            flux_balance_analysis_dict(m, GLPK.Optimizer)["BIOMASS_Ecoli_core_w_GAM"] /
            original_flux["BIOMASS_Ecoli_core_w_GAM"],
)

# (You may notice that restricting the ATP maintenance pseudo-reaction (`ATPM`)
# had a mildly surprising effect of actually increasing the biomass production
# by a few percent.  That is because the cells are not required to produce ATP
# to survive and may invest the nutrients and energy elsewhere.)

# ## Screening with reaction combinations

# The same analysis can be scaled up to screen for combinations of critical
# reactions, giving possibly more insight into the redundancies in the model:

running_reaction_combinations = [
    (rid1, rid2, x1, x2) for (rid1, x1) in running_reactions,
    (rid2, x2) in running_reactions
]

biomass_mtx = screen(
    model,
    variants = [
        [
            with_changed_bound(
                rid1,
                lower_bound = -0.5 * abs(x1),
                upper_bound = 0.5 * abs(x1),
            ),
            with_changed_bound(
                rid2,
                lower_bound = -0.5 * abs(x2),
                upper_bound = 0.5 * abs(x2),
            ),
        ] for (rid1, rid2, x1, x2) in running_reaction_combinations
    ],
    analysis = m ->
        flux_balance_analysis_dict(m, GLPK.Optimizer)["BIOMASS_Ecoli_core_w_GAM"] /
        original_flux["BIOMASS_Ecoli_core_w_GAM"],
)

# Finally, let's plot the result:

using CairoMakie, Clustering

order =
    hclust([
        sum((i .- j) .^ 2) for i in eachcol(biomass_mtx), j in eachcol(biomass_mtx)
    ]).order

labels = first.(running_reactions)[order];
positions = collect(eachindex(labels))

f = Figure(fontsize = 8)
ax = Axis(f[1, 1], xticks = (positions, labels), yticks = (positions, labels))
heatmap!(ax, positions, positions, biomass_mtx[order, order])
ax.xticklabelrotation = π / 3
ax.xticklabelalign = (:right, :center)
ax.yticklabelrotation = π / 6
ax.yticklabelalign = (:right, :center)
f

# Remember that [`screen`](@ref) can be parallelized just [by supplying worker
# IDs](../distributed/1_functions.md). Use that to gain significant speedup
# with analyses of larger models.
