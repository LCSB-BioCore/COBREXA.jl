# # Exploring model contents

# For practical reasons, COBREXA.jl supports many different model types. These
# comprise ones that reflect the storage formats (such as [`JSONModel`](@ref)
# and [`SBMLModel`](@ref)), and ones that are more easily accessible for users
# and mimic the usual workflows in COBRA methodology:
#
# - [`ObjectModel`](@ref), which contains and object-oriented representation
#   of model internals, built out of [`Reaction`](@ref), [`Metabolite`](@ref)
#   and [`Gene`](@ref) structures, in a way similar to e.g.
#   [COBRApy](https://github.com/opencobra/cobrapy/)
# - [`MatrixModel`](@ref), which contains array-oriented representation of the
#   model structures, such as stoichiometry matrix and the bounds vector, in a
#   way similar to e.g. [COBRA
#   toolbox](https://github.com/opencobra/cobratoolbox)

# The fields in [`ObjectModel`](@ref) structure can be discovered using `fieldnames` as follows:

using COBREXA

fieldnames(ObjectModel)

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

sm = load_model(ObjectModel, "e_coli_core.json")
typeof(sm.reactions)

fieldnames(Reaction)

# This process (along with e.g. Tab completion in REPL) allows you to pick
# various information about many objects, for example about a specific
# reaction:

sm.reactions["TALA"].name
#
sm.reactions["TALA"].grr #gene-reaction relationship
#
sm.reactions["TALA"].subsystem
#
sm.reactions["TALA"].ub #upper rate bound

# The same applies to [`MatrixModel`](@ref):

fieldnames(MatrixModel)
#
cm = load_model(MatrixModel, "e_coli_core.json")
#
cm.S
#
cm.rxns[1:10]
