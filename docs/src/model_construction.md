# Model Construction

## Defining genes
Genes are represented by the `Gene` type in `COBREXA.jl`, see [Model
Structure](@ref) for details.  `Gene`s can be constructed using either an empty
constructor, or a constructor taking only the string `id` of the gene.

```@docs
Gene()
Gene(::String)
```

```@example
using COBREXA

gene = Gene("gene1")
gene.name = "gene 1 name"
gene # pretty printing
```

Helper functions from Base have also been overwritten to make accessing arrays
of genes easy.

```@docs
findfirst(::Array{Gene, 1}, ::String)
getindex(::Array{Gene, 1}, ::Gene)
```

## Defining metabolites

Metabolites are represented by the `Metabolite` type in `COBREXA.jl`, see
[Model Structure](@ref) for details.  The simplest way to define a new
metabolite is by using the empty constructor `Metabolite()`.

Alternatively, `Metabolite(id::String)` can be used to assign only the `id`
field of the `Metabolite`.

```@docs
Metabolite()
Metabolite(::String)
```

The other fields can be modified as usual, if desired.

```@example
using COBREXA # hide
atp = Metabolite("atp")
atp.name = "Adenosine triphosphate"
atp.formula = "C10H12N5O13P3"
atp.charge = -4 # note that the charge is an Int

atp # pretty printing
```

Basic analysis of `Metabolite`s is also possible. The molecular structure of a
metabolite can be inspected by calling `get_atoms(met::Metabolite)`.  This
function is useful for checking atom balances across reactions, or the entire
model.

```@docs
get_atoms(met::Metabolite)
```

```@example
using COBREXA # hide
atp = Metabolite("atp") # hide
atp.name = "Adenosine triphosphate" # hide
atp.formula = "C10H12N5O13P3" # hide
get_atoms(atp)
```

Helper functions from Base have also been overwritten to make accessing arrays
of metabolites easy.

```@docs
findfirst(mets::Array{Metabolite, 1}, metid::String)
getindex(mets::Array{Metabolite, 1}, met::Metabolite)
```

## Defining reactions

Reactions are represented by the `Reaction` type in `COBREXA.jl`, see [Model
Structure](@ref) for details.  The simplest way to define a new reaction is by
using the empty constructor `Reaction()`.  All the other fields still need to
be assigned.

```@docs
Reaction()
```

Another option is to use
`Reaction(id::String, metabolites::Dict{Metabolite, Float64}, dir="bidir")`,
which assigns the reaction `id`, the reaction stoichiometry (through the
metabolite dictionary argument), and the directionality of the reaction.  The
remaining fields still need to be assigned, if desired.

```@docs
Reaction(id::String, metabolites::Dict{Metabolite, Float64}, dir="bidir")
```

```@example
using COBREXA # hide
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

See the discussion in [Model Structure](@ref) about how to assign `grr` to a
reaction.

Yet another way of defining a reaction is through overloading of the operators:
`*, +, ∅, ⟶, →, ←, ⟵, ↔, ⟷`.  The longer and shorter arrows mean the same
thing, i.e. `⟶` is the same as `→`, etc.  The other fields of the reaction
still need to be set directly.

```@example
using COBREXA # hide
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

When building exchange, demand, or sinks reactions the `∅` empty metabolite
should be used to indicate that a metabolite is being created or destroyed.

```@example
using COBREXA # hide
adp = Metabolite("adp") # hide
adp.formula = "C10H12N5O10P2" # hide

ex_rxn = ∅ ⟷ adp # exchange reaction
```

It is also possible to check if a reaction is mass balanced by using
`is_mass_balanced(rxn::Reaction)`.  Note, this function requires that all the
metabolites in the reaction have formulas assigned to them to work properly.

```@docs
is_mass_balanced
```

```@example
using COBREXA # hide
atp = Metabolite("atp") # hide
atp.name = "Adenosine triphosphate" # hide
atp.formula = "C10H12N5O13P3" # hide
atp.charge = -4 # hide
adp = Metabolite("adp") # hide
adp.formula = "C10H12N5O10P2" # hide

unbal_rxn = 2.0adp ⟶ 1.0*atp # unbalanced reaction
is_mass_balanced(unbal_rxn)
```

Helper functions from Base have also been overwritten to make accessing arrays
of reactions easy.

```@docs
findfirst(rxns::Array{Reaction, 1}, rxnid::String)
getindex(rxns::Array{Reaction, 1}, rxn::Reaction)
```

## Building models

Once you have defined some metabolites, genes, and reactions, you can construct
a model! This is most simply done by using the empty model constructor:

```@docs
CobraModel()
```

The fields of `CobraModel` can then be assigned as usual.

```@example
using COBREXA

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

anabolism = 10.0 * atp + 10.0*h2o ⟶ 10.0*adp + 10.0*pp
anabolism.id = "anabolism"

glc_ex = ∅ ⟷ glc # exchanges are defined like this so negative fluxes mean import
lac_ex = ∅ ⟷ lac # positive flux means export

mets = [atp, adp, pp, h2o, glc, lac]
genes = [g1, g2, g3, g4]
rxns = [catabolism, anabolism, lac_ex, glc_ex]

model = CobraModel()
model.id = "Test model"
add!(model, mets)
add!(model, rxns)
add!(model, genes)

model # pretty printing
```

Here the `add` functions were used to add the reactions, metabolites and genes
to the model.

```@docs
add!(model::CobraModel, rxns::Array{Reaction, 1})
add!(model::CobraModel, mets::Array{Metabolite, 1})
add!(model::CobraModel, genes::Array{Gene, 1})
```

Checking for duplicates of genes, metabolites or reactions can also be done.
Note that duplicate checking is NOT performed when models are imported.  If you
suspect the model quality you should check each reaction, metabolite and gene
yourself.

```@docs
check_duplicate_annotations(genes::Array{Gene, 1}, gene::Gene)
check_duplicate_annotations(mets::Array{Metabolite, 1}, cmet::Metabolite)
check_duplicate_annotations(rxns::Array{Reaction, 1}, crxn::Reaction)
check_same_formula(mets::Array{Metabolite, 1}, met::Metabolite)
check_duplicate_reaction(rxns::Array{Reaction, 1}, crxn::Reaction)
```

```@example duplex
using COBREXA

met1 = Metabolite()
met1.id = "met1"
met1.name = "Metabolite 1"
met1.formula = "C6H12O6N"
met1.charge = 1
met1.compartment = "c"
met1.notes = Dict("notes"=>["This is a made up metabolite", "Another note"])
met1.annotation = Dict("sboterm" => "sbo000001", "kegg.compound" => ["C0001", "C0010"])

met2 = Metabolite("met2")
met2.formula = "C6H12O6N"

met3 = Metabolite("met3")
met3.formula = "X"
met3.annotation = Dict("sboterm" => "sbo00001", "kegg.compound" => ["C02222", "C0001"])

mets = [met1, met2, met3]

dup, ind = check_duplicate_annotations(mets, met3)
```

```@example duplex
mms = check_same_formula([met3, met1], met2)
```

Similar functionality exists for genes and reactions. Duplicate reactions,
metabolites or genes can be removed using `rm!`.

```@docs
rm!(model::CobraModel, rxns::Union{Array{Reaction, 1}, Reaction})
rm!(model::CobraModel, mets::Union{Array{Metabolite, 1}, Metabolite})
rm!(model::CobraModel, genes::Union{Array{Gene, 1}, Gene})
```

```@example
using COBREXA # hide
atp = Metabolite("atp")
atp2 = Metabolite("atp2")

mets = [atp, atp2]

model = CobraModel()
add!(model, mets)

rm!(model, atp2)
```

A model can also be checked to ensure that no metabolites or genes are missing
relative to the reactions.  `fix_model` also ensures that no extra metabolites
are present.

```@docs
fix_model!(model::CobraModel)
```

```@example
using COBREXA # hide
atp = Metabolite("atp")
adp = Metabolite("adp")

anabolism = 10.0 * atp ⟶ 10.0*adp
anabolism.id = "anabolism"

mets = [atp]
rxns = [anabolism]

model = CobraModel()
model.id = "Test model"
add!(model, mets) # missing adp
add!(model, rxns)

fix_model!(model) # adp added

model # now has 2 metabolites
```

Helper functions from Base have also been overwritten to make accessing
reactions, metabolites and genes easy from a model.

```@docs
getindex(model::CobraModel, rxn::Reaction)
getindex(model::CobraModel, rxn::Metabolite)
getindex(model::CobraModel, rxn::Gene)
```
