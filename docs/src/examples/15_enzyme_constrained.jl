# # GECKO

# GECKO algorithm can be used to easily adjust the metabolic activity within the
# cell to respect many known parameters, measured by proteomics and other
# methods.
#
# The original description from GECKO is by: [Sánchez, et. al., "Improving the
# phenotype predictions of a yeast genome‐scale metabolic model by incorporating
# enzymatic constraints.", Molecular systems biology,
# 2017](https://doi.org/10.15252/msb.20167411).
#
# The analysis method and implementation in COBREXA is similar to
# [sMOMENT](14_simplified_enzyme_constrained.md), but GECKO is able to process and represent much
# larger scale of the constraints -- mainly, it supports multiple isozymes for
# each reaction, and the isozymes can be grouped into "enzyme mass groups" to
# simplify interpretation of data from proteomics.

# For demonstration, we will generate artificial random data in a way similar
# to the [sMOMENT example](14_simplified_enzyme_constrained.md):

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

using COBREXA, GLPK

model = load_model("e_coli_core.json")

import Random
Random.seed!(1) # repeatability

gene_product_masses = Dict(genes(model) .=> randn(n_genes(model)) .* 10 .+ 60)

rxns = filter(
    x ->
        !looks_like_biomass_reaction(x) &&
            !looks_like_exchange_reaction(x) &&
            !isnothing(reaction_gene_associations(model, x)),
    variables(model),
)

# The main difference from sMOMENT comes from allowing multiple isozymes per
# reaction (reactions with missing isozyme informations will be ignored,
# leaving them as-is):
rxn_isozymes = Dict(
    rxn => [
        Isozyme(
            Dict(isozyme_genes .=> 1),
            randn() * 100 + 600, #forward kcat
            randn() * 100 + 500, #reverse kcat
        ) for isozyme_genes in reaction_gene_associations(model, rxn)
    ] for rxn in rxns
)

# We also construct similar bounds for total gene product amounts:
gene_product_bounds = Dict(genes(model) .=> Ref((0.0, 10.0)))

# With this, the construction of the model constrained by all enzymatic
# information is straightforward:

enzyme_constrained_model =
    model |> with_enzyme_constraints(;
        reaction_isozymes = rxn_isozymes,
        gene_product_bounds,
        gene_product_molar_mass = gene_product_masses,
        gene_product_mass_group = _ -> "uncategorized", # all products belong to the same "uncategorized" category
        gene_product_mass_group_bound = _ -> 100.0, # the total limit of mass in the single category
    )

# (Alternatively, you may use [`make_enzyme_constrained_model`](@ref), which does the same
# without piping by `|>`.)

# The stoichiometry and coupling in the enzyme_constrained model is noticeably more complex;
# you may notice new "reactions" added that simulate the gene product
# utilization:

[stoichiometry(enzyme_constrained_model); coupling(enzyme_constrained_model)]

# Again, the resulting model can be used in any type of analysis. For example, flux balance analysis:

opt_model = flux_balance_analysis(enzyme_constrained_model, GLPK.Optimizer)

# Get the fluxes

flux_sol = flux_dict(enzyme_constrained_model, opt_model)

# Get the gene product concentrations

gp_concs = gene_product_dict(enzyme_constrained_model, opt_model)

# Get the total masses assigned to each mass group

gene_product_mass_group_dict(enzyme_constrained_model, opt_model)

# Variability:

flux_variability_analysis(
    enzyme_constrained_model,
    GLPK.Optimizer,
    bounds = gamma_bounds(0.95),
)

# ...and sampling:
affine_hit_and_run(
    enzyme_constrained_model,
    warmup_from_variability(enzyme_constrained_model, GLPK.Optimizer),
)' * reaction_variables(enzyme_constrained_model)
