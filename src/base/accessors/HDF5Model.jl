
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
    _h5_read_sparse(SparseMat, model.h5["stoichiometry"])
end

function bounds(model::HDF5Model)::Tuple{Vector{Float64},Vector{Float64}}
    precache!(model)
    (HDF5.readmmap(model.h5["lower_bounds"]), HDF5.readmmap(model.h5["upper_bounds"]))
end

function balance(model::HDF5Model)::SparseVec
    precache!(model)
    _h5_read_sparse(SparseVec, model.h5["balance"])
end

function objective(model::HDF5Model)::SparseVec
    precache!(model)
    _h5_read_sparse(SparseVec, model.h5["objective"])
end
