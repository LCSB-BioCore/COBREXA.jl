# # sMOMENT

# sMOMENT algorithm can be used to easily adjust the metabolic activity within
# the cell to respect known enzymatic parameters and enzyme mass constraints
# measured by proteomics and other methods.
#
# The original description from sMOMENT is by [Bekiaris, and Klamt, "Automatic
# construction of metabolic models with enzyme constraints.", BMC
# bioinformatics, 2020](https://doi.org/10.1186/s12859-019-3329-9)
#
# Let's load some packages:

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

using COBREXA, GLPK

model = load_model("e_coli_core.json")

# We will necessarily need the enzyme turnover numbers (aka "kcats") and masses
# of the required gene products.  You do not necessarily need to know all data
# for the given model, but the more you have, the better the approximation will
# be.
#
# For the demonstration purpose, we will generate the data randomly. In a
# realistic setting, you would input experimental or database-originating data
# here:

import Random
Random.seed!(1) # repeatability

gene_product_masses = Dict(genes(model) .=> randn(n_genes(model)) .* 10 .+ 60)

# We only take the reactions that have gene products (i.e., enzymes) associated with them):
rxns = filter(
    x ->
        !looks_like_biomass_reaction(x) &&
            !looks_like_exchange_reaction(x) &&
            !isnothing(reaction_gene_association(model, x)),
    reactions(model),
)

# The information about each enzyme and its capabilities is stored in an
# [`Isozyme`](@ref) structure. For simplicity, sMOMENT ignores much of the
# information about the multiplicity of required gene products and
# other.

rxn_isozymes = Dict(
    rxn => Isozyme(
        Dict(vcat(reaction_gene_association(model, rxn)...) .=> 1),
        randn() * 100 + 600, #forward kcat
        randn() * 100 + 500, #reverse kcat
    ) for rxn in rxns
)

# In case some of the reactions are missing in `rxn_isozymes`, sMOMENT simply
# ignores them.
#
# Once the data is gathered, we create a model that wraps the original model
# with additional sMOMENT structure:

smoment_model =
    model |> with_smoment(
        reaction_isozyme = rxn_isozymes,
        gene_product_molar_mass = gene_product_masses,
        total_enzyme_capacity = 50.0,
    )

# (You could alternatively use the [`make_smoment_model`](@ref) to create the
# structure more manually; but [`with_smoment`](@ref) is easier to use e.g.
# with [`screen`](@ref).)

# In turn, you should have a complete model with unidirectional reactions and
# additional coupling, as specified by the sMOMENT method:

[stoichiometry(smoment_model); coupling(smoment_model)]

# the type (SMomentModel) is a model wrapper -- it is a thin additional layer
# that just adds the necessary sMOMENT-relevant information atop the original
# model, which is unmodified. That makes the process very efficient and
# suitable for large-scale data processing. You can still access the original
# "base" model hidden in the SMomentModel using [`unwrap_model`](@ref).

# Other than that, the [`SMomentModel`](@ref) is a model type like any other,
# and you can run any analysis you want on it, such as FBA:

flux_balance_analysis_dict(smoment_model, GLPK.Optimizer)

# (Notice that the total reaction fluxes are reported despite the fact that
# reactions are indeed split in the model! The underlying mechanism is provided
# by [`reaction_flux`](@ref) accessor.)

# [Variability](06_fva.md) of the sMOMENT model can be explored as such:

flux_variability_analysis(smoment_model, GLPK.Optimizer, bounds = gamma_bounds(0.95))

# ...and a sMOMENT model sample can be obtained [as usual with
# sampling](16_hit_and_run.md):

(
    affine_hit_and_run(
        smoment_model,
        warmup_from_variability(smoment_model, GLPK.Optimizer),
    )' * reaction_flux(smoment_model)
)
