
# # Loading models

# `COBREXA` can load models stored in `.mat`, `.json`, and `.xml` formats (with
# the latter denoting SBML formatted models).
#
# We will primarily use the *E. coli* "core" model to demonstrate the utilities
# found in `COBREXA`. First, let's download the model in several formats.

## Downloads the model files if they don't already exist
!isfile("e_coli_core.mat") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.mat", "e_coli_core.mat");
!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");
!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml");

#md # !!! tip "Save bandwidth!"
#md #     The published models usually do not change very often. It is
#md #     therefore pretty useful to save them to a central location and load
#md #     them from there. That saves your time, and does not unnecessarily
#md #     consume the connectivity resources of the model repository.

# Load the models using the [`load_model`](@ref) function. Models are able to
# "pretty-print" themselves, hiding the inner complexity:

using COBREXA

mat_model = load_model("e_coli_core.mat")
#

json_model = load_model("e_coli_core.json")
#

sbml_model = load_model("e_coli_core.xml")
#

#md # !!! note "Note: `load_model` infers the input type from the file extension"
#md #       Notice how each model was read into memory as a model type corresponding
#md #       to its file type, i.e. the file ending with `.json` loaded as a
#md #       [`JSONModel`](@ref), the file ending with `.mat` loaded as [`MATModel`](@ref), and the
#md #       file ending with `.xml` loaded as an [`SBMLModel`](@ref).

# The loaded models contain the data in a format that is preferably as
# compatible as possible with the original representation. In particular, the
# JSON model contains the representation of the JSON tree:

json_model.json

# SBML models contain a complicated structure from [`SBML.jl`
# package](https://github.com/LCSB-BioCore/SBML.jl):

typeof(sbml_model.sbml)

# MAT models contain MATLAB data:

mat_model.mat

# In all cases, you can access the data in the model in the same way, e.g.,
# using [`reactions`](@ref) to get a list of the reactions in the models:

reactions(mat_model)[1:5]
#

reactions(json_model)[1:5]

# You can use the [generic accessors](03_exploring.md) to gather more information about
# the model contents, [convert the models](02_convert_save.md) into formats more suitable for
# hands-on processing, and export them back to disk after the
# modification.
#
# All model types can be directly [used in analysis functions](05a_fba.md), such as
# [`flux_balance_analysis`](@ref).
