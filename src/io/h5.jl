
"""
$(TYPEDSIGNATURES)

Return a HDF5Model associated with the given file. Does not actually load
anything (for efficiency) -- use [`precache!`](@ref) to start pulling data into
the memory.
"""
function load_h5_model(file_name::String)::HDF5Model
    return HDF5Model(file_name)
end

"""
$(TYPEDSIGNATURES)

Converts and writes a metabolic model to disk in the HDF5 format.

Additionally returns an (uncached) [`HDF5Model`](@ref) that represents the
contents of the saved file. Because all HDF5-based models need to be backed by
disk storage, writing the data to disk (using this function) is the only way to
make new HDF5 models.
"""
function save_h5_model(model::MetabolicModel, file_name::String)::HDF5Model
    rxns = reactions(model)
    rxnp = sortperm(rxns)
    mets = metabolites(model)
    metp = sortperm(mets)
    h5open(file_name, "w") do f
        write(f, "metabolites", mets[metp])
        write(f, "reactions", rxns[rxnp])
        _h5_write_sparse(create_group(f, "balance"), balance(model)[metp])
        _h5_write_sparse(create_group(f, "objective"), objective(model)[rxnp])
        _h5_write_sparse(create_group(f, "stoichiometry"), stoichiometry(model)[metp, rxnp])
        let (lbs, ubs) = bounds(model)
            write(f, "lower_bounds", lbs[rxnp])
            write(f, "upper_bounds", ubs[rxnp])
        end
    end
    # TODO: genes, grrs, compartments. Perhaps chemistry and others?
    HDF5Model(file_name)
end

"""
$(TYPEDSIGNATURES)

Close (and un-cache) the [`HDF5Model`](@ref) data. This allows the associated
file to be opened for writing again.
"""
function Base.close(model::HDF5Model)
    if !isnothing(model.h5)
        close(model.h5)
        model.h5 = nothing
    end
end


_h5_mmap_nonempty(x) = length(x) > 0 ? HDF5.readmmap(x) : HDF5.read(x)

function _h5_write_sparse(g::HDF5.Group, v::SparseVector)
    write(g, "n", v.n)
    write(g, "nzind", v.nzind)
    write(g, "nzval", v.nzval)
end

function _h5_read_sparse(::Type{X}, g::HDF5.Group) where {X<:SparseVector}
    n = read(g["n"])
    nzind = _h5_mmap_nonempty(g["nzind"])
    nzval = _h5_mmap_nonempty(g["nzval"])
    SparseVector{eltype(nzval),eltype(nzind)}(n, nzind, nzval)
end

function _h5_write_sparse(g::HDF5.Group, m::SparseMatrixCSC)
    write(g, "m", m.m)
    write(g, "n", m.n)
    write(g, "colptr", m.colptr)
    write(g, "rowval", m.rowval)
    write(g, "nzval", m.nzval)
end

function _h5_read_sparse(::Type{X}, g::HDF5.Group) where {X<:SparseMatrixCSC}
    m = read(g["m"])
    n = read(g["n"])
    colptr = _h5_mmap_nonempty(g["colptr"])
    rowval = _h5_mmap_nonempty(g["rowval"])
    nzval = _h5_mmap_nonempty(g["nzval"])
    SparseMatrixCSC{eltype(nzval),eltype(colptr)}(m, n, colptr, rowval, nzval)
end

