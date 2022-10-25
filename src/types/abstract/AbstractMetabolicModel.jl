
"""
    abstract type AbstractMetabolicModel end

A helper supertype of everything usable as a linear-like model for COBREXA
functions.

If you want your model type to work with COBREXA, add the `AbstractMetabolicModel` as
its supertype, and implement the accessor functions. Accessors
[`reactions`](@ref), [`metabolites`](@ref), [`stoichiometry`](@ref),
[`bounds`](@ref) and [`objective`](@ref) must be implemented; others are not
mandatory and default to safe "empty" values.
"""
abstract type AbstractMetabolicModel end

"""
    abstract type ModelWrapper <: AbstractMetabolicModel end

A helper supertype of all "wrapper" types that contain precisely one other
[`AbstractMetabolicModel`](@ref).
"""
abstract type ModelWrapper <: AbstractMetabolicModel end

const SparseMat = SparseMatrixCSC{Float64,Int}
const SparseVec = SparseVector{Float64,Int}
const MatType = AbstractMatrix{Float64}
const VecType = AbstractVector{Float64}
const StringVecType = AbstractVector{String}

"""
    GeneAssociation = Vector{Vector{String}}

An association to genes, represented as a logical formula in a positive
disjunctive normal form (DNF). (The 2nd-level vectors of strings are connected
by "and" to form conjunctions, and the 1st-level vectors of these conjunctions
are connected by "or" to form the DNF.)
"""
const GeneAssociation = Vector{Vector{String}}

"""
    MetaboliteFormula = Dict{String,Int}

Dictionary of atoms and their abundances in a molecule.
"""
const MetaboliteFormula = Dict{String,Int}

"""
    Annotations = Dict{String,Vector{String}}

Dictionary used to store (possible multiple) standardized annotations of
something, such as a [`Metabolite`](@ref) and a [`Reaction`](@ref).

# Example
```
Annotations("PubChem" => ["CID12345", "CID54321"])
```
"""
const Annotations = Dict{String,Vector{String}}

"""
    Notes = Dict{String,Vector{String}}

Free-form notes about something (e.g. a [`Gene`](@ref)), categorized by
"topic".
"""
const Notes = Dict{String,Vector{String}}
