
# # Converting, modifying and saving models

# COBREXA.jl can export JSON and MATLAB-style model formats, which can be
# useful when exchanging the model data with other software.
#
# For a test, let's download and open a SBML model:

using COBREXA

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml");

sbml_model = load_model("e_coli_core.xml")


# You can save the model as `.json` or `.mat` file using the
# [`save_model`](@ref) function:

save_model(sbml_model, "converted_e_coli.json")
save_model(sbml_model, "converted_e_coli.mat")


# ## Using serialization for quick loading and saving

# If you are saving the models only for future processing in Julia environment,
# it is often wasteful to encode the models to external formats and decode them
# back. Instead, you can use the "native" Julia data format, accessible with
# package `Serialization`.
#
# This way, you can use `serialize` to save any model format (even the
# complicated [`ObjectModel`](@ref), which does not have a "native" file format
# representation):

using Serialization

sm = convert(ObjectModel, sbml_model)

open(f -> serialize(f, sm), "myModel.stdmodel", "w")

# The models can then be loaded back using `deserialize`:

sm2 = deserialize("myModel.stdmodel")
issetequal(metabolite_ids(sm), metabolite_ids(sm2))

# This form of loading operation is usually pretty quick:
t = @elapsed deserialize("myModel.stdmodel")
@info "Deserialization took $t seconds"
# Notably, large and complicated models with thousands of reactions and
# annotations can take tens of seconds to decode properly. Serialization allows
# you to minimize this overhead, and scales well to tens of millions of
# reactions.

#md # !!! warning "Compatibility"
#md #     The format of serialized models may change between Julia versions.
#md #     In particular, never use the the serialized format for publishing models -- others will have hard time finding the correct Julia version to open them.
#md #     Similarly, never use serialized models for long-term storage -- your future self will have hard time finding the historic Julia version that was used to write the data.

# ## Converting and saving a modified model

# To modify the models easily, it is useful to convert them to a format that
# simplifies this modification. You may use e.g. [`MatrixModel`](@ref) that
# exposes the usual matrix-and-vectors structure of models as used in MATLAB
# COBRA implementations, and [`ObjectModel`](@ref) that contains structures,
# lists and dictionaries of model contents, as typical in Python COBRA
# implementations. The object-oriented nature of [`ObjectModel`](@ref) is
# better for making small modifications that utilize known identifiers of model
# contents.
#
# Conversion of any model to [`ObjectModel`](@ref) can be performed using the
# standard Julia `convert`:

sm = convert(ObjectModel, sbml_model)

# The conversion can be also achieved right away when loading the model, using
# an extra parameter of [`load_model`](@ref):

sm = load_model(ObjectModel, "e_coli_core.json")

# As an example, we change an upper bound on one of the reactions:

sm.reactions["PFK"].ub = 10.0

# After [possibly applying more modifications](04_standardmodel.md), you can again save the
# modified model in a desirable exchange format:

save_model(sm, "modified_e_coli.json")
save_model(sm, "modified_e_coli.mat")

# More information about [`ObjectModel`](@ref) internals is available [in a
# separate example](04_standardmodel.md).
