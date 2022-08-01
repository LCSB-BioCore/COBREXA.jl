# # Loopless FBA

# Here we will use [`flux_balance_analysis`](@ref) and
# [`flux_variability_analysis`](@ref) to analyze a toy model of *E. coli* that
# is constrained in a way that removes all thermodynamically infeasible loops in the flux solution.
# For more details about the algorithm, see Schellenberger, Lewis, and, Palsson. "Elimination
# of thermodynamically infeasible loops in steady-state metabolic models.", *Biophysical
# journal*, 2011 (https://doi.org/10.1016/j.bpj.2010.12.3707).

# If it is not already present, download the model:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model("e_coli_core.xml")

# In COBREXA.jl, the Loopless FBA is implemented as a modification of the
# normal FBA, called [`add_loopless_constraints`](@ref).
# We use GLPK optimizer here, because the loopless constraints add integer
# programming into the problem. Simpler solvers (such as Tulip) may not be able
# to solve the mixed integer-linear (MILP) programs.

loopless_flux = flux_balance_analysis_vec(
    model,
    GLPK.Optimizer,
    modifications = [add_loopless_constraints()],
)

# The representation is particularly convenient since it allows to also explore
# other properties of loopless models, such as variability and parsimonious
# balance, as well as other analyses that accept `modifications` parameter:

loopless_variability = flux_variability_analysis(
    model,
    GLPK.Optimizer,
    modifications = [add_loopless_constraints()],
)

# For details about the loopless method, refer to Schellenberger, Jan, Nathan
# E. Lewis, and Bernhard Ø. Palsson: "Elimination of thermodynamically
# infeasible loops in steady-state metabolic models." *Biophysical journal*
# 100, no. 3 (2011), pp. 544-553.
