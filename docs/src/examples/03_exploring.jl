# # Exploring model contents

# For practical reasons, COBREXA.jl supports many different model types. These comprise ones that reflect the storage formats (such as [`JSONModel`](@ref) and [`SBMLModel`](@ref)), and ones that are more easily accessible for users and mimic the usual workflows in COBRA methodology:
#
# - [`StandardModel`](@ref), which contains and object-oriented representation
# of model internals, built out of [`Reaction`](@ref), [`Metabolite`](@ref) and
# [`Gene`](@ref) structures, in a way similar to e.g.
# [COBRApy](https://github.com/opencobra/cobrapy/)
# - [`CoreModel`](@ref), which contains array-oriented representation of the
# model structures, such as stoichiometry matrix and the bounds vector, in a
# way similar to e.g. [COBRA
# toolbox](https://github.com/opencobra/cobratoolbox)

# The fields in [`StandardModel`](@ref) structure can be discovered using `fieldnames` as follows:

fieldnames(StandardModel)

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

sm = load_model(StandardModel, "e_coli_core.json")
typeof(sm.reactions)

fieldnames(Reaction)

# This process (along with e.g. Tab completion in REPL) allows you to pick
# various information about many objects, for example about a specific
# reaction:

sm.reactions["TALA"].name

sm.reactions["TALA"].grr #gene-reaction relationship

sm.reactions["TALA"].subsystem

sm.reactions["TALA"].ub #upper rate bound

# The same applies to [`CoreModel`](@ref):

fieldnames(CoreModel)

cm = load_model(CoreModel, "e_coli_core.json")

cm.S

cm.rxns[1:10]

# ## Generic accessors

# To prevent the complexities of object representation, `COBREXA.jl` uses a set
# of generic interface functions that can extract various important information
# from any supported model type. This approach ensures that the analysis
# functions can work on any data.

# For example, you can check the reactions and metabolites contained in SBML
# and JSON models using the same accessor:

reactions(sm)

reactions(cm)

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

[x.name for x in methodswith(MetabolicModel, COBREXA) if endswith(String(x.file), "MetabolicModel.jl")]
