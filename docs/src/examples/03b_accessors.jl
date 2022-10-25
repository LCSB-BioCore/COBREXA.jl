
# # Generic accessors

# To prevent the complexities of object representation, `COBREXA.jl` uses a set
# of generic interface functions that can extract various important information
# from any supported model type. This approach ensures that the analysis
# functions can work on any data.

# For example, you can check the reactions and metabolites contained in any
# model type ([`SBMLModel`](@ref), [`JSONModel`](@ref), [`CoreModel`](@ref),
# [`StandardModel`](@ref), and any other) using the same accessor:

using COBREXA

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

js = load_model("e_coli_core.json")
reactions(js)
#
std = convert(CoreModel, js)
reactions(std)

# All accessors allow systematic access to information about reactions,
# stoichiometry, metabolite properties and chemistry, genes, and various model
# annotations.
#
# The most notable ones include:
#
# - [`reactions`](@ref), [`metabolites`](@ref) and [`genes`](@ref) return
#   respective vectors of identifiers of reactions, metabolites and genes present
#   in the model,
# - [`stoichiometry`](@ref) returns the S matrix
# - [`balance`](@ref) returns the right-hand vector of the linear model in form `Ax=b`
# - [`bounds`](@ref) return lower and upper bounds of reaction rates
# - [`metabolite_charge`](@ref) and [`metabolite_formula`](@ref) return details about metabolites
# - [`objective`](@ref) returns the objective of the model (usually labeled as `c`)
# - [`reaction_gene_association`](@ref) describes the dependency of a reaction on gene products
#
# A complete, up-to-date list of accessors can be always generated using `methodswith`:

using InteractiveUtils

accessors = [
    x.name for x in methodswith(AbstractMetabolicModel, COBREXA) if
    endswith(String(x.file), "AbstractMetabolicModel.jl")
]

println.(accessors);

#md # !!! note "Note: Not all accessors may be implemented for all the models"
#md #       It is possible that not all the accessors are implemented for all the model
#md #       types. If this is the case, usually `nothing` or an empty data structure is
#md #       returned. If you need a specific accessor, just overload the function you require!
