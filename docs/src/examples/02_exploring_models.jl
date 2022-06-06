# # Accessing model internals using the generic accessors

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# ## Using the generic interface to access model details

# To prevent the complexities of object representation, `COBREXA.jl` uses a set
# of generic interface functions that extract various important information
# from all supported model types. This approach ensures that the analysis
# functions can work on any data.

# For example, you can check the reactions and metabolites contained in SBML
# and JSON models using the same accessor:

reactions(json_model)
#

reactions(sbml_model)
#

issetequal(reactions(json_model), reactions(mat_model)) # do models contain the same reactions?

# All accessors are defined in a single file in COBREXA source code; you may
# therefore get a list of all accessors as follows:
using InteractiveUtils

for method in filter(
    x -> endswith(string(x.file), "MetabolicModel.jl"),
    InteractiveUtils.methodswith(MetabolicModel, COBREXA),
)
    println(method.name)
end