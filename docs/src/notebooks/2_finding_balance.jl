# # Finding balance and variability of constraint-based models

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

using Tulip, OSQP

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

# ## Modifications

# Often it is desirable to add a slight modififaction to the problem before
# performing analysis, to see e.g. differences of the model behavior caused by
# the change introduced.
#
# `COBREXA.jl` supports several modifications by default, which include
# changing objective sense, optimizer attributes, flux constraints,
# optimization objective, reaction and gene knockouts, and others.

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
        silence, # and make it quiet.
    ],
)

# ## Flux balance analysis with molecular crowding (FBAwMC)

# Flux balance with molecular crowding incorporates enzyme capacity constraints into the
# classic flux balance analysis algorithm. Essentially, an extra constraint is added to
# the optimization problem: `∑ wᵢ × vᵢ ≤ 1` where `wᵢ` weights each internal flux `vᵢ`. See
# `Beg, Qasim K., et al. "Intracellular crowding defines the mode and sequence of
# substrate uptake by Escherichia coli and constrains its metabolic activity." Proceedings of
# the National Academy of Sciences 104.31 (2007): 12663-12668.` for more details.

# First load the model
model = load_model("e_coli_core.json")

# Next, simulate the model over a range of substrate uptake rates.
without_crowding = Dict{Float64,Vector{Float64}}()
with_crowding = Dict{Float64,Vector{Float64}}()
glucose_uptakes = collect(-(1.0:0.5:20))

for glc in glucose_uptakes
    no_crowding = flux_balance_analysis_dict( # classic FBA
        model,
        Tulip.Optimizer;
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
            change_constraint("EX_glc__D_e"; lb = glc),
        ],
    )

    without_crowding[glc] =
        [no_crowding["BIOMASS_Ecoli_core_w_GAM"], no_crowding["EX_ac_e"]]

    crowding = flux_balance_analysis_dict( # FBAwMC
        model,
        Tulip.Optimizer;
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 1000),
            add_crowding_constraint(0.004), # crowding constraint gets added here
            change_constraint("EX_glc__D_e"; lb = glc),
        ],
    )

    with_crowding[glc] = [crowding["BIOMASS_Ecoli_core_w_GAM"], crowding["EX_ac_e"]]
end

# Finally, plot the results to compare classic FBA with FBAwMC.
using CairoMakie
fig = Figure();
ax1 = Axis(fig[1, 1]);
lines!(
    ax1,
    -glucose_uptakes,
    [without_crowding[glc][1] for glc in glucose_uptakes],
    label = "no crowding",
    linewidth = 5,
    linestyle = :dash,
)
lines!(
    ax1,
    -glucose_uptakes,
    [with_crowding[glc][1] for glc in glucose_uptakes],
    label = "with crowding",
    linewidth = 5,
    linestyle = :dot,
)
ax1.ylabel = "Growth rate [1/h]"
ax2 = Axis(fig[2, 1])
lines!(
    ax2,
    -glucose_uptakes,
    [without_crowding[glc][2] for glc in glucose_uptakes],
    label = "no crowding",
    linewidth = 5,
    linestyle = :dash,
)
lines!(
    ax2,
    -glucose_uptakes,
    [with_crowding[glc][2] for glc in glucose_uptakes],
    label = "with crowding",
    linewidth = 5,
    linestyle = :dot,
)
fig[1:2, 2] = Legend(fig, ax1, "Constraint")
ax2.xlabel = "Glucose uptake rate [mmol/gDW/h]"
ax2.ylabel = "Acetate production rate\n[mmol/gDW/h]"
fig

# Notice that overflow metabolism is modeled by incorporating the crowding constraint.

#md # !!! tip "Tip: code your own modification like [`add_crowding`](@ref)"
#md #       Making custom problem modification functions is really simple due to the
#md #       tight intergration between COBREXA and JuMP. Look at the source code for
#md #       the implemented modifications in `src\analysis\modifications` to get a flavour
#md #       for it.
