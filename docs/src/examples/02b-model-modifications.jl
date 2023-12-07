
# # Making adjustments to the model
#
# Typically, we do not need to solve the models as they come from the authors
# (someone else already did that!), but we want to perform various
# perturbations in the model structure and conditions, and explore how the
# model behaves in the changed conditions.
#
# With COBREXA, there are 2 different approaches that one can take:
# 1. We can change the model structure and use the changed metabolic model.
# This is better for doing simple and small but systematic modifications, such
# as removing metabolites, adding reactions, etc.
# 2. We can intercept the pipeline that converts the metabolic model to
# constraints and then to the optimizer representation, and make small
# modifications along that way. This is better for various technical model
# adjustments, such as using combined objectives or adding reaction-coupling
# constraints.
#
# Here we demonstrate the first, "modelling" approach. The main advantage of
# that approach is that the modified model is still a FBC model, and you can
# export, save and share it via the AbstractFBCModels interace. The main
# disadvantage is that the "common" FBC model interface does not easily express
# various complicated constructions (communities, reaction coupling, enzyme
# constraints, etc.) -- see the [example about modifying the
# constraints](02c-constraint-modifications.md) for a closer look on how to
# modify even such complex constructions.
#
# TODO here. :)
