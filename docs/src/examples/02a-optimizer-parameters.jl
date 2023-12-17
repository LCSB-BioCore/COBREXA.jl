
# # Changing optimizer parameters
#
# Many optimizers require fine-tuning to produce best results. You can pass in
# additional optimizer settings via the `modifications` parameter of
# [`flux_balance`](@ref). These include e.g.
#
# - [`set_optimizer_attribute`](@ref) (typically allowing you to tune e.g.
#   iteration limits, tolerances, or floating-point precision)
# - [`set_objective_sense`](@ref) (allowing you to change and reverse the
#   optimization direction, if required)
# - [`silence`](@ref) to disable the debug output of the optimizer
# - and even [`set_optimizer`](@ref), which changes the optimizer
#   implementation used (this is not quite useful in this case, but becomes
#   beneficial with more complex, multi-stage optimization problems)
#
# To demonstrate this, we'll use the usual toy model:

using COBREXA
import JSONFBCModels, Tulip

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

model = load_model("e_coli_core.json")

# Running a FBA with a silent optimizer that has slightly increased iteration
# limit for IPM algorithm may now look as follows:
solution = flux_balance_analysis(
    model,
    Tulip.Optimizer;
    modifications = [silence, set_optimizer_attribute("IPM_IterationsLimit", 1000)],
)

@test !isnothing(solution) #src

# To see some of the effects of the configuration changes, you may e.g.
# deliberately cripple the optimizer's possibilities to a few iterations, which
# will cause it to fail, return no solution, and verbosely describe what
# happened:

solution = flux_balance_analysis(
    model,
    Tulip.Optimizer;
    modifications = [set_optimizer_attribute("IPM_IterationsLimit", 2)],
)

println(solution)

@test isnothing(solution) #src

# Applicable optimizer attributes are documented in the documentations of the
# respective optimizers. To browse the possibilities, you may want to see the
# [JuMP documentation page that summarizes the references to the available
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).
