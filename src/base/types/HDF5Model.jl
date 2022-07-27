
"""
    struct HDF5Model

A model that is stored in HDF5 format. The model data is never really pulled
into memory, but instead mmap'ed as directly as possible into the Julia
structures.  This makes reading the `HDF5Model`s extremely fast.

TODO: it can be sent around like SerializedModel.

TODO: decide how to create the models
- there should be save_model support
- there could be explicit function like `serialize_model`, like `hdf5size_model`? :D
- explicitly say you can't convert to HDF5model.

TODO: decide if this should be able to hold eeeeeeeeverything or just the
CoreModel-style computation-relevant parts.
"""
mutable struct HDF5Model <: MetabolicModel
    h5::Maybe{HDF5.File}
    filename::String

    HDF5Model(filename::String) = new(nothing, filename)
end

function precache!(model::HDF5Model)::Nothing
    if isnothing(model.h5)
        model.h5 = h5open(model.filename, "r")
    end
    nothing
end

function n_reactions(model::HDF5Model)::Int
    precache!(model)
    length(model.h5["reactions"])
end

function reactions(model::HDF5Model)::Vector{String}
    precache!(model)
    # TODO is there any reasonable method to mmap strings from HDF5?
    read(model.h5["reactions"])
end

function n_metabolites(model::HDF5Model)::Int
    precache!(model)
    length(model.h5["metabolites"])
end

function metabolites(model::HDF5Model)::Vector{String}
    precache!(model)
    read(model.h5["metabolites"])
end

function stoichiometry(model::HDF5Model)::SparseMat
    precache!(model)
    h5_read_sparse(SparseMat, model.h5["stoichiometry"])
end

function bounds(model::HDF5Model)::Tuple{Vector{Float64},Vector{Float64}}
    precache!(model)
    (HDF5.readmmap(model.h5["lower_bounds"]), HDF5.readmmap(model.h5["upper_bounds"]))
end

function balance(model::HDF5Model)::SparseVec
    precache!(model)
    h5_read_sparse(SparseVec, model.h5["balance"])
end

function objective(model::HDF5Model)::SparseVec
    precache!(model)
    h5_read_sparse(SparseVec, model.h5["objective"])
end
