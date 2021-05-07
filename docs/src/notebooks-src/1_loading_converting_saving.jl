# # Loading, converting, and saving models

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# `COBREXA` can load models stored in `.mat`, `.json`, and `.xml` formats (with
# the latter denoting SBML formatted models). 

# We will primarily use the *E. coli* toy model to demonstrate the utilities
# found in `COBREXA`. To begin, first download some models in a variety of
# formats.

## Downloads the model files if they don't already exist
!isfile("e_coli_core.mat") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.mat", "e_coli_core.mat")
!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

#md # !!! tip "Tip: Model files" 
#md #     To avoid having to download the models multiple times, save the models in
#md #     a central location and always load them from that location.

#nb # To avoid having to download the models multiple times, save the models in
#nb # a central location and always load them from that location.

# ## Loading models

# Now load the models using the `load_model` function. Each model has "pretty printing".

mat_model = load_model("e_coli_core.mat")
#
json_model = load_model("e_coli_core.json")
#
sbml_model = load_model("e_coli_core.xml")
#

#md # !!! note "Note: Model type is inferred from the file type"
#md #       Notice how each model was read into memory as a model type corresponding
#md #       to its file type, i.e. the file ending with `.json` loaded as a
#md #       `JSONModel`, the file ending with `.mat` loaded as a `MATModel`, and the
#md #       file ending with `.xml` loaded as an `SBMLModel`. 
#md #   
#md #       When loading a model from file like this (i.e., when the file type matches
#md #       model type), no information gets lost because the model is loaded directly
#md #       into memory from the underlying file type.

#nb # Notice how each model was read into memory as a model type corresponding
#nb # to its file type, i.e. the file ending with `.json` loaded as a
#nb # `JSONModel`, the file ending with `.mat` loaded as a `MATModel`, and the
#nb # file ending with `.xml` loaded as an `SBMLModel`. 

#nb #  When loading a model from file like this (i.e., when the file type matches
#nb #  model type), no information gets lost because the model is loaded directly
#nb #  into memory from the underlying file type.

# You can directly inspect the model object, although this is discouraged because
# generic interface functions, which are applicable to __all__ model types,
# are provided by `COBREXA`.

json_model.json # try to not use this functionality

# ## Using the generic interface to access model details

# `COBREXA` makes use of a generic interface that enforces a uniform set of
# accessor functions for model attributes. This approach ensures that any custom
# function you create and test using, e.g. `JSONModel`s, also works on, e.g.
# `SBMLModel`s, etc.

# Here is a list of all the currently supported accessor functions.
using InteractiveUtils #hide 
for (i, method) in enumerate(
    filter(
        x -> endswith(string(x.file), "MetabolicModel.jl"),
        InteractiveUtils.methodswith(MetabolicModel, COBREXA),
    ),
)
    println(i, ") ", method.name)
end

# By using this interface you can access internal details about models using the
# exact same function, regardless of the underlying model.

reactions(json_model)
#
issetequal(reactions(json_model), reactions(mat_model)) # same reactions returned

# ## Converting between models

# It is also possible to convert model types to-and-fro. To do this, use the
# `convert` function, which is overloaded from Julia's `Base` module and functions
# in the same way. 

#md # !!! danger "Data loss may occur when converting between models"
#md #       The generic interface is used to convert between model types, so only
#md #       data accessible through the generic accessors will be converted.

#nb # The generic interface is used to convert between model types, so only
#nb # data accessible through the generic accessors will be converted.


json_model_to_mat_model = convert(MATModel, json_model)

# ## Loading and converting models at the same time

# `COBREXA` also allows you to load models into a specific format using a
# slightly more advanced version of the `load_model` function.

mat_model_from_json_model_file = load_model(MATModel, "e_coli_core.json") # specify the model type to load

# It is important to remember that this format change is always associated with
# a call to `convert`, thus data loss may occur.

# ## Saving models

# `COBREXA` also allows you to save models. The format the models should be
# saved in is inferred from the file suffix.

save_model(mat_model_from_json_model_file, "converted_model.json")
#
save_model(mat_model_from_json_model_file, "converted_model.mat")
