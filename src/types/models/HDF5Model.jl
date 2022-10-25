
"""
$(TYPEDEF)

A model that is stored in HDF5 format. The model data is never really pulled
into memory, but instead mmap'ed as directly as possible into the Julia
structures.  This makes reading the `HDF5Model`s extremely fast, at the same
time the (uncached) `HDF5Model`s can be sent around efficiently among
distributed nodes just like [`Serialized`](@ref) models, provided the nodes
share a common storage.

All HDF5Models must have the backing disk storage. To create one, use
[`save_h5_model`](@ref) or [`save_model`](@ref) with `.h5` file extension. To
create a temporary model that behaves like a model "in memory", save it to a
temporary file. For related reasons, you can not use `convert` models to
`HDF5Model` format, because the conversion would impliy having the model saved
somewhere.

# Fields
$(TYPEDFIELDS)
"""
mutable struct HDF5Model <: AbstractMetabolicModel
    h5::Maybe{HDF5.File}
    filename::String

    HDF5Model(filename::String) = new(nothing, filename)
end

function Accessors.precache!(model::HDF5Model)::Nothing
    if isnothing(model.h5)
        model.h5 = h5open(model.filename, "r")
    end
    nothing
end

function Accessors.n_reactions(model::HDF5Model)::Int
    precache!(model)
    length(model.h5["reactions"])
end

function Accessors.reactions(model::HDF5Model)::Vector{String}
    precache!(model)
    # TODO is there any reasonable method to mmap strings from HDF5?
    read(model.h5["reactions"])
end

function Accessors.n_metabolites(model::HDF5Model)::Int
    precache!(model)
    length(model.h5["metabolites"])
end

function Accessors.metabolites(model::HDF5Model)::Vector{String}
    precache!(model)
    read(model.h5["metabolites"])
end

function Accessors.stoichiometry(model::HDF5Model)::SparseMat
    precache!(model)
    h5_read_sparse(SparseMat, model.h5["stoichiometry"])
end

function Accessors.bounds(model::HDF5Model)::Tuple{Vector{Float64},Vector{Float64}}
    precache!(model)
    (HDF5.readmmap(model.h5["lower_bounds"]), HDF5.readmmap(model.h5["upper_bounds"]))
end

function Accessors.balance(model::HDF5Model)::SparseVec
    precache!(model)
    h5_read_sparse(SparseVec, model.h5["balance"])
end

function Accessors.objective(model::HDF5Model)::SparseVec
    precache!(model)
    h5_read_sparse(SparseVec, model.h5["objective"])
end
