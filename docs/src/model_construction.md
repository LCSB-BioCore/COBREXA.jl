# Model Construction

## Defining genes
Genes are represented by the `Gene` type in `CobraTools`, see [Model Structure](@ref) for details.
`Gene`s can be constructed using either an empty constructor, or a constructor taking only
the string `id` of the gene.
```@docs
Gene()
Gene(::String)
```
```@example
using CobraTools # hide

gene = Gene("gene1")
gene.name = "gene 1 name"
gene # pretty printing
```
Helper functions from Base have also been overwritten to make accessing arrays of genes easy.
```@@docs
findfirst(::Array{Gene, 1}, ::String)
getindex(::Array{Gene, 1}, ::Gene)
```
## Defining metabolites
Metabolites are represented by the `Metabolite` type in `CobraTools`, see [Model Structure](@ref) for details. 
The simplest way to define a new metabolite is by using the empty constructor `Metabolite()`. 
Alternatively, `Metabolite(id::String)` can be used to assign only the `id` field of the `Metabolite`. 
```@docs
Metabolite()
Metabolite(::String)
```
The other fields can be modified as usual, if desired.
```@example
using CobraTools

atp = Metabolite("atp")
atp.name = "Adenosine triphosphate"
atp.formula = "C10H12N5O13P3"
atp.charge = -4 # note that the charge is an Int

atp # pretty printing
```

Basic analysis of `Metabolite`s is also possible. The molecular structure of a metabolite can be inspected by calling
`get_atoms(met::Metabolite)`.
This function is useful for checking atom balances across reactions, or the entire model.
```@docs
CobraTools.get_atoms(met::Metabolite)
```
```@example
using CobraTools # hide
atp = Metabolite("atp") # hide
atp.name = "Adenosine triphosphate" # hide
atp.formula = "C10H12N5O13P3" # hide

get_atoms(atp)
```
Helper functions from Base have also been overwritten to make accessing arrays of metabolites easy.
```@@docs
findfirst(mets::Array{Metabolite, 1}, metid::String)
getindex(mets::Array{Metabolite, 1}, met::Metabolite)
```
## Defining reactions
Reactions are represented by the `Reaction` type in `CobraTools`, see [Model Structure](@ref) for details.
The simplest way to define a new reaction is by using the empty constructor `Reaction()`. 
All the other fields still need to be assigned.
```@docs
Reaction()
```

Another option is to use `Reaction(id::String, metabolites::Dict{Metabolite, Float64}, dir="bidir")`, which 
assigns the reaction `id`, the reaction stoichiometry (through the metabolite dictionary argument), and the directionality of the reaction.
The remaining fields still need to be assigned, if desired.
```@docs
Reaction(id::String, metabolites::Dict{Metabolite, Float64}, dir="bidir")
```
```@example
using CobraTools # hide
atp = Metabolite("atp") # hide
atp.name = "Adenosine triphosphate" # hide
atp.formula = "C10H12N5O13P3" # hide
atp.charge = -4 # hide
gene = Gene("gene1") # hide

adp = Metabolite("adp") # define another metabolite

metdict = Dict(atp => -1.0, adp => 1.0) # nb stoichiometries need to be floats

rxn = Reaction("dummy rxn", metdict, "for")
rxn.annotation["ec-code"] = ["0.0.0.0"]
rxn.grr = [[gene]] # only gene1 is required for this reaction to work

rxn # pretty printing
```
See the discussion in [Model Structure](@ref) about how to assign `grr` to a reaction. 

Yet another way of defining a reaction is through overloading of the operators: `*, +, ∅, ⟶, →, ←, ⟵, ↔, ⟷`.
The longer and shorter arrows mean the same thing, i.e. `⟶` is the same as `→`, etc.
The other fields of the reaction still need to be set directly.
```@example
using CobraTools # hide
atp = Metabolite("atp") # hide
atp.name = "Adenosine triphosphate" # hide
atp.formula = "C10H12N5O13P3" # hide
atp.charge = -4 # hide
adp = Metabolite("adp") # hide
adp.formula = "C10H12N5O10P2" # hide

another_rxn = 2.0adp ⟶ 2.0*atp # forward reaction
another_rxn.id = "another dummy rxn"
another_rxn
```
When building exchange, demand, or sinks reactions the `∅` empty metabolite should be used to indicate that a metabolite is being created or destroyed.
```@example
using CobraTools # hide
adp = Metabolite("adp") # hide
adp.formula = "C10H12N5O10P2" # hide

ex_rxn = ∅ ⟷ adp # exchange reaction
```
It is also possible to check if a reaction is mass balanced by using `is_mass_balanced(rxn::Reaction)`. 
Note, this function requires that all the metabolites in the reaction have formulas assigned to them to work properly.
```@docs
is_mass_balanced
```
```@example
using CobraTools # hide
atp = Metabolite("atp") # hide
atp.name = "Adenosine triphosphate" # hide
atp.formula = "C10H12N5O13P3" # hide
atp.charge = -4 # hide
adp = Metabolite("adp") # hide
adp.formula = "C10H12N5O10P2" # hide

unbal_rxn = 2.0adp ⟶ 1.0*atp # unbalanced reaction
is_mass_balanced(unbal_rxn)
```
Helper functions from Base have also been overwritten to make accessing arrays of reactions easy.
```@@docs
findfirst(rxns::Array{Reaction, 1}, rxnid::String)
getindex(rxns::Array{Reaction, 1}, rxn::Reaction)
```
## Building models
Once you have defined some metabolites, genes, and reactions, you can construct a model! This is most simply done by
using the empty model constructor:
```@docs
Model()
```
The fields of `CobraTools.Model` can then be assigned as usual.
```@example
using CobraTools

atp = Metabolite("atp")
adp = Metabolite("adp")
glc = Metabolite("glc")
h2o = Metabolite("h2o")
pp = Metabolite("pp") # pi (phosphate) usually means π, so pp is used here instead
lac = Metabolite("lac")

g1 = Gene("gene1")
g2 = Gene("gene2")
g3 = Gene("gene3")
g4 = Gene("gene4")

catabolism = glc + 2.0*adp + 2.0*pp ⟶ 2.0 * lac + 2.0 * h2o + 2.0 * atp
catabolism.id = "catabolism"
catabolism.grr = [[g1, g2], [g3, g4]]

anabolism = 10.0 * atp ⟶ 10.0*adp
anabolism.id = "anabolism"

glc_ex = ∅ ⟷ glc # exchanges are defined like this so negative fluxes mean import
lac_ex = ∅ ⟷ lac # positive flux means export

mets = [atp, adp, pp, h2o, glc, lac]
genes = [g1, g2, g3, g4]
rxns = [catabolism, anabolism, lac_ex, glc_ex]

model = Model()
model.id = "Test model"
add!(model, mets)
add!(model, rxns)
add!(model, genes)

model # pretty printing
```
Here the `add` functions were used to add the reactions, metabolites and genes to the model.
```@docs
add!(model::CobraTools.Model, rxns::Array{Reaction, 1})
add!(model::CobraTools.Model, mets::Array{Metabolite, 1})
add!(model::CobraTools.Model, genes::Array{Gene, 1})
```
Checking for duplicates can also be done. 
Note that duplicate checking is NOT performed when models are imported. 
If you suspect the model quality you should check each reaction, metabolite and gene yourself.
```@docs
is_duplicate(model::CobraTools.Model, rxn::Reaction)
is_duplicate(model::CobraTools.Model, cmet::Metabolite)
is_duplicate(model::CobraTools.Model, cgene::Gene)
```
```@example
using CobraTools # hide
atp = Metabolite("atp") 
adp = Metabolite("adp") 

anabolism = 10.0 * atp ⟶ 10.0*adp
anabolism.id = "anabolism"

anabolism2 = 10.0 * atp ⟶ 10.0*adp
anabolism.id = "anabolism2"

mets = [atp, adp]
rxns = [anabolism]

model = Model()
model.id = "Test model"
add!(model, mets)
add!(model, rxns)

is_duplicate(model, anabolism2)
```
Duplicate reactions, metabolites or genes can be removed using `rm!`.
```@docs
rm!(model::CobraTools.Model, rxns::Union{Array{Reaction, 1}, Reaction})
rm!(model::CobraTools.Model, mets::Union{Array{Metabolite, 1}, Metabolite})
rm!(model::CobraTools.Model, genes::Union{Array{Gene, 1}, Gene})
```
```@example
using CobraTools # hide
atp = Metabolite("atp")
atp2 = Metabolite("atp2") 

mets = [atp, atp2]

model = Model()
add!(model, mets)

rm!(model, atp2)
```
A model can also be checked to ensure that no metabolites or genes are missing relative to the reactions. 
`fix_model` also ensures that no extra metabolites are present.
```@docs
fix_model!(model::CobraTools.Model)
```
```@example
using CobraTools # hide
atp = Metabolite("atp") 
adp = Metabolite("adp") 

anabolism = 10.0 * atp ⟶ 10.0*adp
anabolism.id = "anabolism"

mets = [atp]
rxns = [anabolism]

model = Model()
model.id = "Test model"
add!(model, mets) # missing adp
add!(model, rxns)

fix_model!(model) # adp added
```
Helper functions from Base have also been overwritten to make accessing reactions, metabolites and genes easy from a model.
```@@docs
getindex(model::CobraTools.Model, rxn::Reaction)
getindex(model::CobraTools.Model, rxn::Metabolite)
getindex(model::CobraTools.Model, rxn::Gene)
```