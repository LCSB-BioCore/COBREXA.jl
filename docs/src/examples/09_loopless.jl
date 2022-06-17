# # Loopless FBA

# Here we will use [`flux_balance_analysis`](@ref), [`flux_variability_analysis`](@ref),
# [`parsimonious_flux_balance_analysis`](@ref), and
# [`minimize_metabolic_adjustment_analysis`](@ref), along with the modification functions of
# `COBREXA.jl`, to analyze a toy model of *E. coli*.

# If it is not already present, download the model:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model("e_coli_core.xml")

# In COBREXA.jl, the Loopless FBA is implemented as a modification of the
# normal FBA, called [`add_loopless_constraints`](@ref).

loopless_flux = flux_balance_analysis_vec(model, GLPK.Optimizer,
    modifications = [add_loopless_constraints()])

# The representation is particularly convenient since it allows to also explore
# other properties of loopless models, such as variability and parsimonious
# balance, as well as other analyses that accept `modifications` parameter:

loopless_variability = flux_variability_analysis(model, GLPK.Optimizer,
    modifications = [add_loopless_constraints()])
