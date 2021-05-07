# # Constraint-based analysis of single cell models

# In this tutorial we will use `COBREXA`'s `flux_balance_analysis`,
# `flux_variability_analysis`, and `parsimonious_flux_balance_analysis`
# functions to analyze a toy model of *E. coli*.

# If it is not already present, load the model.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")
#
using COBREXA

#md # !!! tip "Tip: use `?` when unsure about a function" 
#md #       When you are unsure about how a function works, use `?
#md #       function_name` to access the docstrings of a function.

#nb # When you are unsure about how a function works, use `? function_name` to
#nb # access the docstrings of a function.


model = load_model("e_coli_core.xml")

# ## Optimization solvers in `COBREXA`
#
# To actually perform any optimization based analysis we need to load an
# optimizer. Any [`JuMP` supported
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# will work. Here we will use [`Tulip.jl`](https://github.com/ds4dm/Tulip.jl) to
# optimize linear programs and
# [OSQP.jl](https://osqp.org/docs/get_started/julia.html) to optimize quadratic
# programs.

using Tulip
using OSQP

lp_optimizer = Tulip.Optimizer
qp_optimizer = OSQP.Optimizer

# ## Flux balance analysis
#
# Most analysis functions come in two flavours that determine how the output
# will be returned. The suffix of the function determines this behaviour. 
# 1. `flux_balance_analysis_vec` returns a vector of fluxes in the order of the
#       reactions returned by `reactions` (recall that this is one of the generic
#       interface accessors). 
# 2. `flux_balance_analysis_dict` returns a
#       dictionary mapping reaction ids to fluxes. 
# In both cases there are two required inputs: the `model` and the `optimizer`.

vec_soln = flux_balance_analysis_vec(model, lp_optimizer)
#
dict_soln = flux_balance_analysis_dict(model, lp_optimizer)

# ## Problem modification

# Often it is desirable to modify the problem before performing analysis.
# Problem modifications include things like changing the: objective sense,
# optimizer, solver attributes, flux constraints, and optimization objective. It
# is also possible to knock out genes. For completeness we demonstrate all of
# the mentioned problem modifications below.

#md # !!! note "Note: modifications added to analysis functions are temporary"  
#md #       `COBREXA` supports two flavours of modifications to models. The type
#md #       shown here is temporary, i.e. after the analysis function terminates the model
#md #       is the same as it was before. See the reconstruction tutorial for more information
#md #       about ways to permanently change the model.

#nb # Note: `COBREXA` supports two flavours of modifications to models. The type
#nb # shown here is temporary, i.e. after the analysis function terminates the model
#nb # is the same as it was before. See the reconstruction tutorial for more information
#nb # about ways to permanently change the model.

dict_soln = flux_balance_analysis_dict(
    model,
    qp_optimizer; # will change to Tulip below
    modifications = [ # modifications are evaluated in order
        change_objective("R_BIOMASS_Ecoli_core_w_GAM"),
        change_constraint("R_EX_glc__D_e", -12, -12),
        knockout(["b0978", "b0734"]), # knocks out cytochrome oxidase (CYTBD) 
        change_optimizer(lp_optimizer), # swap back to using Tulip
        change_optimizer_attribute("IPM_IterationsLimit", 110), # this is a Tulip specific attribute, other solvers have other attributes
        change_sense(MAX_SENSE), # another valid option is MIN_SENSE
    ],
)

# ## Flux variability analysis

# Flux variability analysis can be performed in serial or in parallel. Here the
# serial version is demonstrated, see the more advanced tutorials for the
# parallel implementation.

fva_mins, fva_maxs = flux_variability_analysis_dict(
    model,
    Tulip.Optimizer;
    bounds = objective_bounds(0.99), # the objective function is allowed to vary by ~1% from the FBA optimum
    modifications = [
        change_solver_attribute("IPM_IterationsLimit", 500),
        change_constraint("R_EX_glc__D_e", -10, -10),
        change_constraint("R_EX_o2_e", 0.0, 0.0),
    ],
)

fva_maxs["R_EX_ac_e"]["R_EX_ac_e"] # acetate exchange maximized, acetate exchange flux

# ## Parsimonious flux balance analysis

# Parsimonious flux balance analysis finds a unique flux solution that minimizes
# the sum of fluxes of the system subject to maintaining the same objective
# value as the flux balance analysis solution. This function requires the use of
# a quadratic program optimizer (`OSQP.jl` is used in this example). Like in
# `flux_balance_analysis`, two variants exist where the suffix determines the
# function output.

#md # !!! note  "A quadratic programming optimizer is required" 
#md #           If you are using an optimizer that supports quadratic programming, like
#md #           `Gurobi.jl`, then you only need to specify one optimizer e.g. 
#md #           `parsimonious_flux_balance_analysis_dict(model, Gurobi.Optimizer)`.
#md #           All problem and/or optimizer modifications should be passed to the 
#md #           `modifications` keyword argument. The `qp_modifications` keyword
#md #           argument should be ignored in this case. If you are using two optimizers,
#md #           then the optimizer passed first is used as the linear
#md #           programming optimizer, and the second as the quadratic
#md #           programming optimizer. This is demonstrated below.

#nb # If you are using an optimizer that supports quadratic programming, like
#nb # `Gurobi.jl`, then you only need to specify one optimizer e.g. 
#nb # `parsimonious_flux_balance_analysis_dict(model, Gurobi.Optimizer)`.
#nb #  All problem and/or optimizer modifications should be passed to the 
#nb # `modifications` keyword argument. The `qp_modifications` keyword
#nb # argument should be ignored in this case. If you are using two optimizers,
#nb # then the optimizer passed first is used as the linear
#nb # programming optimizer, and the second as the quadratic
#nb # programming optimizer. This is demonstrated below.

dict_soln = parsimonious_flux_balance_analysis_dict(
    model,
    lp_optimizer;
    modifications = [
        change_constraint("R_EX_glc__D_e", -12, -12),
        change_solver_attribute("IPM_IterationsLimit", 500), # modifies the required optimizer input (lp_optimizer in this case)
    ],
    qp_modifications = [
        change_solver(qp_optimizer), # only necessary if the first solver cannot handle QPs
        change_solver_attribute("verbose", false),
    ],
)
#
vec_soln = parsimonious_flux_balance_analysis_vec(
    model,
    lp_optimizer;
    modifications = [
        change_constraint("R_EX_glc__D_e", -12, -12),
        change_solver_attribute("IPM_IterationsLimit", 500), # modifies the required optimizer input
    ],
    qp_modifications = [
        change_solver(qp_optimizer), # only necessary if the first solver cannot handle QPs
        change_solver_attribute("verbose", false),
    ],
)
