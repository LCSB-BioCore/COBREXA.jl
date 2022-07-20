# # Parsimonious flux balance analysis (pFBA)

# Parsimonious flux balance analysis attempts to find a realistic flux of a
# model, by trying to minimize squared sum of all fluxes while maintaining the
# reached optimum. COBREXA.jl implements it in function
# [`parsimonious_flux_balance_analysis`](@ref) (accompanied by vector- and
# dictionary-returning variants
# [`parsimonious_flux_balance_analysis_vec`](@ref) and
# [`parsimonious_flux_balance_analysis_dict`](@ref)).
#
# As usual, we demonstrate the functionality on the *E. Coli* model:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, Tulip, OSQP

model = load_model("e_coli_core.xml")

# Because the parsimonious objective is quadratic, we need a proper optimizer.
# As the simplest choice, we can use
# [`OSQP.jl`](https://osqp.org/docs/get_started/julia.html), but any any
# [`JuMP.jl`-supported
# optimizer](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# that supports quadratic programming will work.

#md # !!! note "Note: OSQP can be sensitive"
#md #       We recommend reading the documentation of `OSQP` before using it, since
#md #       it may give inconsistent results depending on what settings
#md #       you use. Commercial solvers like `Gurobi`, `Mosek`, `CPLEX`, etc.
#md #       require less user engagement.

# Running of basic pFBA is perfectly analogous to running of [FBA](05a_fba.md)
# and other analyses. We add several modifications that improve the solution
# (using functions [`silence`](@ref), and
# [`change_optimizer_attribute`](@ref)), and fix the glucose exchange (using
# [`change_constraint`](@ref)) in order to get a more reasonable result:

fluxes = parsimonious_flux_balance_analysis_dict(
    model,
    OSQP.Optimizer;
    modifications = [
        silence, # optionally silence the optimizer (OSQP is very verbose by default)
        change_optimizer_attribute("polish", true), # tell OSQP to invest time into improving the precision of the solution
        change_constraint("R_EX_glc__D_e"; lb = -12, ub = -12), # fix glucose consumption rate
    ],
)

# ## Using different optimizers for linear and quadratic problems
#
# It is quite useful to use specialized optimizers for specialized tasks in
# pFBA. In particular, one would usually require to get a precise solution from
# the linear programming part (where the precision is achievable), and trade
# off a little precision for vast improvements in computation time in the
# quadratic programming part.
#
# In pFBA, we can use the `modifications` and `qp_modifications` parameters to
# switch and parametrize the solvers in the middle of the process, which allows
# us to implement precisely that improvement. We demonstrate the switching on a
# vector-returning variant of pFBA:

flux_vector = parsimonious_flux_balance_analysis_vec(
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
