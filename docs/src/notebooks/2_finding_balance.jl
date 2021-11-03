# # Finding balance and variability of constraint-based models

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we use [`flux_balance_analysis`](@ref),
# [`flux_variability_analysis`](@ref), and
# [`parsimonious_flux_balance_analysis`](@ref) of `COBREXA.jl` functions to
# analyze a toy model of *E. coli*.

# If it is not already present, download the model.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

#md # !!! tip "Tip: use `?` to get quick help about functions"
#md #       When you are unsure about how a function works, write `?
#md #       function_name` to see the function reference documentation.

model = load_model("e_coli_core.xml")

# ## Optimization solvers in `COBREXA`
#
# To actually perform any optimization based analysis we need to load an
# optimizer. Any [`JuMP.jl`-supported
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# will work. Here, we will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl)
# to optimize linear programs and
# [`OSQP.jl`](https://osqp.org/docs/get_started/julia.html) to optimize quadratic
# programs.

#md # !!! note "Note: OSQP can be sensitive"
#md #       We recommend reading the docs of `OSQP` before using it, since
#md #       it may give inconsistent results depending on what settings
#md #       you use. Commercial solvers like `Gurobi`, `Mosek`, `CPLEX`, etc.
#md #       require less user engagement.

using Tulip, OSQP, GLPK

# ## Flux balance analysis (FBA)
#
# Most analysis functions come in several variants that produce different types
# of output. All of them usually require a model and `JuMP.jl`-compatible
# optimizer to work in the model.
#
# In the case of FBA, you may choose from these variants (here using the
# `Tulip` optimizer):

vec_soln = flux_balance_analysis_vec(model, Tulip.Optimizer)
#
dict_soln = flux_balance_analysis_dict(model, Tulip.Optimizer)

# ## Extended FBA through modifications

# Often it is desirable to add a slight modification to the problem before
# performing analysis, to see e.g. differences of the model behavior caused by
# the change introduced.
#
# `COBREXA.jl` supports several modifications by default, which include changing objective
# sense, optimizer attributes, flux constraints, optimization objective, reaction and gene
# knockouts, and others. These modifications are applied in the order they are specified. It
# is up to the user to ensure that the changes are sensible.

dict_soln = flux_balance_analysis_dict(
    model,
    OSQP.Optimizer;
    modifications = [ # modifications are applied in order
        ## this changes the objective to maximize the biomass production
        change_objective("R_BIOMASS_Ecoli_core_w_GAM"),

        ## this fixes a specific rate of the glucose exchange
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),

        ## this knocks out two genes, i.e. constrains their associated reactions to zero.
        knockout(["b0978", "b0734"]), ## the gene IDs are cytochrome oxidase (CYTBD)

        ## ignore the optimizer specified above and change it to Tulip
        change_optimizer(Tulip.Optimizer),

        ## set a custom attribute of the Tulip optimizer (see Tulip docs for more possibilities)
        change_optimizer_attribute("IPM_IterationsLimit", 110),

        ## explicitly tell the optimizer to maximize the new objective
        change_sense(MAX_SENSE),
    ],
)

# This solution can be display using `flux_summary`. Note, this pretty printing only works
# on flux solutions that are represented as dictionaries.
flux_summary(dict_soln)

# ## Flux variability analysis (FVA)

# The default FVA in [`flux_variability_analysis`](@ref) returns maximized and
# minimized reaction fluxes in a matrix. Here we use the dictionary variant in
# flux_variability_analysis_dict, to show how to easily access specific fluxes
# from its results.

fva_mins, fva_maxs = flux_variability_analysis_dict(
    model,
    Tulip.Optimizer;
    bounds = objective_bounds(0.99), # the objective function is allowed to vary by ~1% from the FBA optimum
    modifications = [
        change_optimizer_attribute("IPM_IterationsLimit", 500),
        change_constraint("R_EX_glc__D_e"; lb = -10, ub = -10),
        change_constraint("R_EX_o2_e"; lb = 0.0, ub = 0.0),
    ],
)

#
fva_maxs["R_EX_ac_e"]["R_EX_ac_e"] # get the maximal acetate exchange flux

# Another option is to display this information using `flux_variability_summary`. This
# pretty printing only works on flux variability analysis results where dictionary keys indicate
# which flux is optimized and the associated value is a flux dictionary.
flux_variability_summary((fva_mins, fva_maxs))

# More sophisticated variants of [`flux_variability_analysis`](@ref) can be used to extract
# specific pieces of information from the solved optimization problems. Here the objective
# value of the minimized flux and the associated biomass growth rate is returned instead
# of every flux.

biomass_idx = first(indexin(["R_BIOMASS_Ecoli_core_w_GAM"], reactions(model))) # index of biomass function
vs = flux_variability_analysis(
    model,
    Tulip.Optimizer;
    bounds = objective_bounds(0.50), # biomass can vary up to 50% less than optimum
    modifications = [
        change_optimizer_attribute("IPM_IterationsLimit", 500),
        change_constraint("R_EX_glc__D_e"; lb = -10, ub = -10),
        change_constraint("R_EX_o2_e"; lb = 0.0, ub = 0.0),
    ],
    ret = m ->
        (COBREXA.JuMP.objective_value(m), COBREXA.JuMP.value(m[:x][biomass_idx])), # m is the model and m[:x] extracts the fluxes from the model
)
#
fva_mins = Dict(rxn => flux for (rxn, flux) in zip(reactions(model), vs[:, 1]))

# ## Parsimonious flux balance analysis (pFBA)

# Parsimonious flux balance analysis (here in
# [`parsimonious_flux_balance_analysis`](@ref) finds a unique flux solution
# that minimizes the squared sum of fluxes of the system subject, while
# maintaining the same objective value as the flux balance analysis solution.
# Since we are optimizing a quadratic objective, we also need to switch to a
# quadratic optimizer. In this case, OSQP will work. We demonstrate it on the
# dictionary-returning variant of pFBA,
# [`parsimonious_flux_balance_analysis_dict`](@ref):

dict_soln = parsimonious_flux_balance_analysis_dict(
    model,
    OSQP.Optimizer;
    modifications = [
        silence, # silence the optimizer (OSQP is very verbose by default)
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),
    ],
)

# The function also has the expectable second variant that returns a vector of
# solutions, in [`parsimonious_flux_balance_analysis_vec`](@ref). Here, we
# utilize it to show how to use different optimizers for finding the optimum
# and for solving the quadratic problem. That may be preferable if the
# optimizer qualities differ for the differing tasks. pFBA allows you to
# specify `qp_modifications` that are applied after the original optimum is
# found, and before the quadratic part of the problem solving begins.

vec_soln = parsimonious_flux_balance_analysis_vec(
    model,
    Tulip.Optimizer; # start with Tulip
    modifications = [
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12),
        change_optimizer_attribute("IPM_IterationsLimit", 500), # we may change Tulip-specific attributes here
    ],
    qp_modifications = [
        change_optimizer(OSQP.Optimizer), # now switch to OSQP (Tulip wouldn't be able to finish the computation)
        change_optimizer_attribute("polish", true), # get an accurate solution, see OSQP's documentation
        silence, # and make it quiet.
    ],
)

# ## Minimizing metabolic adjustment analysis (MOMA)

# MOMA is a technique used to find a flux distribution that is closest to some reference
# distribution with respect to the Euclidian norm.

reference_fluxes = parsimonious_flux_balance_analysis_dict( # reference distribution
    model,
    OSQP.Optimizer;
    modifications = [silence, change_optimizer_attribute("polish", true)],
)

moma = minimize_metabolic_adjustment_analysis_dict(
    model,
    reference_fluxes,
    OSQP.Optimizer;
    modifications = [
        silence,
        change_optimizer_attribute("polish", true),
        change_constraint("R_CYTBD"; lb = 0.0, ub = 0.0), # find flux distribution closest to the CYTBD knockout
    ],
)

# ## Composing (more complicated) modifications

# `COBREXA.jl` contains a number of modifications that allow the user to analyze
# non-standard variants of the classic FBA problem. These include thermodynamic
# ([`add_loopless_constraints`](@ref)), capacity ([`add_crowding_constraints`](@ref)), and
# kinetic/capacity ([`add_moment_constraints`](@ref)) modifications. The documentation of
# each modification details what their purpose is. Here, we will demonstrate how these
# modifications can be composed to generate even more interesting analyses.

# Download the json formatted model so that the reaction identifiers correspond to the Escher map.
!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

model = load_model("e_coli_core.json")

# First find a flux distribution that is thermodyamically loopless and incorporates enzyme
# capacity constraints by composing loopless FBA and FBAwMC.

rid_crowding_weight = Dict(# crowding needs a weight for each flux
    rid => 0.004 for rid in reactions(model) if
    !looks_like_biomass_reaction(rid) && !looks_like_exchange_reaction(rid)
)

loopless_crowding_fluxes = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer;
    modifications = [
        add_crowding_constraints(rid_crowding_weight),
        add_loopless_constraints(),
    ],
)

# Next find a flux distribution that satisfies kinetic/capacity constraints using the moment
# algorithm that is closest (using MOMA) to the loopless crowding solution.

ksas = Dict(rid => 1000.0 for rid in reactions(model)) # make up specific activities of the enyzmes
protein_mass_fraction = 0.56

moment_moma = minimize_metabolic_adjustment_analysis_dict(
    model,
    loopless_crowding_fluxes,
    OSQP.Optimizer;
    modifications = [
        silence,
        change_optimizer_attribute("polish", true),
        change_constraint("EX_glc__D_e", lb = -1000),
        change_constraint("CYTBD"; lb = 0, ub = 0),
        add_moment_constraints(ksas, protein_mass_fraction;),
    ],
)

# Finally, plot the results to inspect the flux distributions visually
using CairoMakie, Escher, ColorSchemes

!isfile("e_coli_core_map.json") && download(
    "http://bigg.ucsd.edu/escher_map_json/e_coli_core.Core%20metabolism",
    "e_coli_core_map.json",
)


maxflux = maximum(abs.(values(moment_moma)))
minflux = minimum(abs.(values(moment_moma)))

# Scale color of reaction edges to fluxes (manually binned)
color_interp(x) = begin
    normed_x = (abs(x) - minflux) / (maxflux - minflux)
    if 0 <= normed_x < 0.01
        ColorSchemes.RdYlBu_4[4]
    elseif 0.01 <= normed_x < 0.25
        ColorSchemes.RdYlBu_4[3]
    elseif 0.25 <= normed_x < 0.5
        ColorSchemes.RdYlBu_4[2]
    else
        ColorSchemes.RdYlBu_4[1]
    end
end
rc = Dict(k => color_interp(v) for (k, v) in moment_moma) # map reaction id to reaction edge color

fig = Figure();
ax = Axis(fig[1, 1]);
escherplot!(ax, "e_coli_core_map.json", reaction_edge_colors = rc)
hidexdecorations!(ax)
hideydecorations!(ax)
fig

#md # !!! tip "Tip: code your own modification like [`add_crowding_constraints`](@ref)"
#md #       Making custom problem modification functions is really simple due to the
#md #       tight intergration between COBREXA and JuMP. Look at the source code for
#md #       the implemented modifications in `src\analysis\modifications` to get a flavour
#md #       for it.
