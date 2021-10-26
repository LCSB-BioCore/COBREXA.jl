
# # Using a custom model data structure

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# This notebooks shows how to utilize the generic accessors and modification
# functions in COBREXA.jl to run the analysis on any custom model type. We will
# create a simple dictionary-style structure that describes the model, allow
# COBREXA to run a FVA on it, and create a simple reaction-removing
# modification.

# First, let's define a very simple stoichiometry-only structure for the model:

using COBREXA

mutable struct MyReaction
    max_rate::Float64 # maximum absolute conversion rate
    stoi::Dict{String,Float64} # stoichimetry of the reaction

    MyReaction() = new(0.0, Dict{String,Float64}())
end

mutable struct MyModel <: MetabolicModel
    optimization_target::String # the "objective" reaction name
    reactions::Dict{String,MyReaction} # dictionary of reactions

    MyModel() = new("", Dict{String,MyReaction}())
    MyModel(o, r) = new(o, r)
end

# With this, we can start defining the accessors:

COBREXA.n_reactions(m::MyModel) = length(m.reactions)
COBREXA.reactions(m::MyModel) = sort(collect(keys(m.reactions)))

# Metabolites are defined only very implicitly, so let's just make a function
# that collects all names. [`n_metabolites`](@ref) can be left at the default
# definition that just measures the output of [`metabolites`](@ref).

function COBREXA.metabolites(m::MyModel)
    mets = Set{String}()
    for (_, r) in m.reactions
        for (m, _) in r.stoi
            push!(mets, m)
        end
    end
    return sort(collect(mets))
end

# Now, the extraction of the linear model. Remember the order of element in the
# vectors must match the order in the output of [`reactions`](@ref) and
# [`metabolites`](@ref).

using SparseArrays

function COBREXA.bounds(m::MyModel)
    max_rates = [m.reactions[r].max_rate for r in reactions(m)]
    (sparse(-max_rates), sparse(max_rates))
end

function COBREXA.objective(m::MyModel)
    if m.optimization_target in keys(m.reactions)
        c = spzeros(n_reactions(m))
        c[first(indexin([m.optimization_target], reactions(m)))] = 1.0
        c
    else
        throw(
            DomainError(
                m.optimization_target,
                "The target reaction for flux optimization not found",
            ),
        )
    end
end

function COBREXA.stoichiometry(m::MyModel)
    sparse([
        get(m.reactions[rxn].stoi, met, 0.0) for met in metabolites(m), rxn in reactions(m)
    ])
end

# Now the model is complete! We can generate a random one and see how it
# performs
import Random
Random.seed!(123)

rxn_names = ["Reaction $i" for i = 'A':'Z'];
metabolite_names = ["Metabolite $i" for i = 1:20];

m = MyModel();
for i in rxn_names
    m.reactions[i] = MyReaction()
end

for i = 1:50
    rxn = rand(rxn_names)
    met = rand(metabolite_names)
    m.reactions[rxn].stoi[met] = rand([-3, -2, -1, 1, 2, 3])
    m.reactions[rxn].max_rate = rand()
end

# Let's see what the model looks like now:
m

# We can run most of the standard function on the model data right away:
using Tulip
m.optimization_target = "Reaction A"
flux_balance_analysis_dict(m, Tulip.Optimizer)


# To be able to use the model conveniently in functions such as
# [`screen`](@ref), you usually want to be able to easily specify the
# modifications. In this example, we enable use of
# [`with_removed_reactions`](@ref) by overloading the internal
# [`remove_reactions`](@ref) for this specific model type:
#
# We need to make an as-shallow-as-possible copy of the model that allows us to
# remove the reactions without breaking the original model.

function COBREXA.remove_reactions(m::MyModel, rxns::AbstractVector{String})
    m = MyModel(m.optimization_target, copy(m.reactions))
    delete!.(Ref(m.reactions), rxns)
    return m
end

# The screening is ready now!

reactions_to_remove = ("Reaction $i" for i = 'B':'Z')

reactions_to_remove .=> screen_variants(
    m,
    [[with_removed_reactions([r])] for r in reactions_to_remove],
    m -> flux_balance_analysis_dict(m, Tulip.Optimizer)["Reaction A"],
)
