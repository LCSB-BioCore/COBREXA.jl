# # Flux balance analysis

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we will use [`flux_balance_analysis`](@ref), [`flux_variability_analysis`](@ref),
# [`parsimonious_flux_balance_analysis`](@ref), and
# [`minimize_metabolic_adjustment_analysis`](@ref), along with the modification functions of
# `COBREXA.jl`, to analyze a toy model of *E. coli*.

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

# ## Extending FBA with modifications

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
