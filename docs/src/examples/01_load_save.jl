# # Loading, converting, and saving models

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# `COBREXA` can load models stored in `.mat`, `.json`, and `.xml` formats (with
# the latter denoting SBML formatted models).

# We will primarily use the *E. Coli* "core" model to demonstrate the utilities
# found in `COBREXA`. First, let's download the model in several formats.

## Downloads the model files if they don't already exist
!isfile("e_coli_core.mat") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.mat", "e_coli_core.mat");
!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");
!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml");

# Now, load the package:

using COBREXA

#md # !!! tip "Save bandwidth!"
#md #     The published models usually do not change very often. It is
#md #     therefore pretty useful to save them to a central location and load
#md #     them from there. That saves your time, and does not unnecessarily
#md #     consume the connectivity resources of the model repository.

# Load the models using the [`load_model`](@ref) function. Each model is able to
# "pretty-print" itself, hiding the inner complexity.

mat_model = load_model("e_coli_core.mat")
#

json_model = load_model("e_coli_core.json")
#

sbml_model = load_model("e_coli_core.xml")
#

#md # !!! note "Note: `load_model` infers the output type from the file extension"
#md #       Notice how each model was read into memory as a model type corresponding
#md #       to its file type, i.e. the file ending with `.json` loaded as a
#md #       [`JSONModel`](@ref), the file ending with `.mat` loaded as [`MATModel`](@ref), and the
#md #       file ending with `.xml` loaded as an [`SBMLModel`](@ref).

# You can directly inspect the model objects, although only with a specific way
# for each specific type.

# JSON models contain their corresponding JSON:

json_model.json

# SBML models contain a complicated structure from [`SBML.jl`
# package](https://github.com/LCSB-BioCore/SBML.jl):

typeof(sbml_model.sbml)

# MAT models contain MATLAB data:

mat_model.mat

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

# ## Converting between model types

# It is possible to convert model types to-and-fro. To do this, use the
# `convert` function, which is overloaded from Julia's `Base`.

#md # !!! danger "Data loss may occur when converting between models"
#md #     The conversion of models only uses the data accessible through the
#md #     generic accessors. Other data may get lost.

m = convert(MATModel, json_model)

# `m` will now contain the MATLAB-style matrix representation of the model:

Matrix(m.mat["S"])

# The loading and conversion can be combined using a shortcut:

m = load_model(MATModel, "e_coli_core.json")

# ## Saving and exporting models

# `COBREXA.jl` supports exporting the models in JSON and MAT format, using [`save_model`](@ref).

save_model(m, "converted_model.json")
save_model(m, "converted_model.mat")

# If you need a non-standard suffix, use the type-specific saving functions:

save_json_model(m, "file.without.a.good.suffix")
save_mat_model(m, "another.file.matlab")

# If you are saving the models only for future processing in Julia environment,
# it is often wasteful to encode the models to external formats and decode them
# back. Instead, you can use the "native" Julia data format, accessible with
# package `Serialization`.
#
# This way, you can use `serialize` to save even the [`StandardModel`](@ref)
# that has no file format associated:

using Serialization

sm = convert(StandardModel, m)

open(f -> serialize(f, sm), "myModel.stdmodel", "w")

# The models can then be loaded back using `deserialize`:

sm2 = deserialize("myModel.stdmodel")
issetequal(metabolites(sm), metabolites(sm2))

# This form of loading operation is usually pretty quick:
t = @elapsed deserialize("myModel.stdmodel")
@info "Deserialization took $t seconds"
# Notably, large and complicated models with thousands of reactions and
# annotations can take seconds to decode properly. Serialization allows you to
# almost completely remove this overhead, and scales well to tens of millions
# of reactions.
