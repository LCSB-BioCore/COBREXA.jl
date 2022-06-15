# # Gene knockouts

# Here we will use the [`knockout`](@ref) function to modify the optimization
# model before solving, in order to simulate genes knocked out. You can pass
# [`knockout`](@ref) to many analysis functions that support parameter
# `modifications`, including [`flux_balance_analysis`](@ref),
# [`flux_variability_analysis`](@ref), and others.

# ## Deleting a single gene

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model("e_coli_core.xml")

# First, let's compute the "original" flux, with no knockouts.
original_flux = flux_balance_analysis_dict(model, GLPK.Optimizer);

# You can find gene IDs that you can knock out using [`genes`](@ref) and
# [`gene_name`](@ref) functions:
genes(model)
# It is possible to sort the genes by gene name to allow easier lookups:
sort(gene_name.(Ref(model), genes(model)) .=> genes(model))

# Compute the flux with a genes knocked out:
flux_with_knockout =
    flux_balance_analysis_dict(model, GLPK.Optimizer, modifications = [knockout("G_b3236")])

# We can see there is a small decrease in production upon knocking out the gene:
biomass_id = "R_BIOMASS_Ecoli_core_w_GAM"
flux_with_knockout[biomass_id] / original_flux[biomass_id]

# Similarly, you can explore how the flux variability has changed once the gene
# is knocked out:
variability_with_knockout =
    flux_variability_analysis(model, GLPK.Optimizer, modifications = [knockout("G_b3236")])

# ## Knocking out multiple genes

# Multiple genes can be knocked out by simply passing a vector of genes to the
# knockout modification. This knocks out all genes that can run the FBA
# reaction:

reaction_gene_association(model, "R_FBA")
#
flux_with_double_knockout = flux_balance_analysis_dict(
    model,
    GLPK.Optimizer,
    modifications = [knockout(["G_b2097", "G_b1773", "G_b2925"])],
)
#
flux_with_double_knockout[biomass_id] / original_flux[biomass_id]

# ## Processing all single gene knockouts
#
# Function [`screen`](@ref) provides a parallelizable and extensible way to run
# the flux balance analysis with the knockout over all genes:

knockout_fluxes = screen(
    model,
    args = tuple.(genes(model)),
    analysis = (m, gene) -> begin
        res = flux_balance_analysis_dict(m, GLPK.Optimizer, modifications = [knockout(gene)])
        if !isnothing(res)
            res[biomass_id]
        end
    end,
)

# It is useful to display the biomass growth rates of the knockout models
# together with the gene name:
sort(gene_name.(Ref(model), genes(model)) .=> knockout_fluxes, by = first)

# ## Processing all multiple-gene deletions
#
# ### Double gene knockouts
#
# Since you can generate any kind of argument matrix for [`screen`](@ref) to
# process, it is straightforward to generate the matrix of all double gene
# knockouts and let the function process it. This computes the biomass
# production of all double-gene knockouts:

gene_groups = [[g1, g2] for g1 in genes(model), g2 in genes(model)];
double_knockout_fluxes = screen(
    model,
    args = tuple.(gene_pairs),
    analysis = (m, gene_groups) -> begin
        res = flux_balance_analysis_dict(
            m,
            GLPK.Optimizer,
            modifications = [knockout(gene_groups)],
        )
        if !isnothing(res)
            res[biomass_id]
        end
    end,
)

# The results can be converted to an easily scrutinizable form as follows:
reshape([gene_name.(Ref(model), p) for p in gene_pairs] .=> double_knockout_fluxes, :)

# ### Triple gene knockouts (and others)
#
# You can extend the same analysis to triple or other gene knockouts by
# generating a different array of gene pairs. For example, you can generate
# gene_groups for triple gene deletion screening:
gene_groups = [[g1, g2, g3] for g1 in genes(model), g2 in genes(model), g3 in genes(model)];

# !!! warning Full triple gene deletion analysis may take a long time to compute.
#     You may use parallel processing with [`screen`](@ref) to speed up the
#     analysis. Alternatively, process only a subset of the genes triples.
