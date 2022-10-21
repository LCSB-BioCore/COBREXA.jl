
h5_mmap_nonempty(x) = length(x) > 0 ? HDF5.readmmap(x) : HDF5.read(x)

function h5_write_sparse(g::HDF5.Group, v::SparseVector)
    write(g, "n", v.n)
    write(g, "nzind", v.nzind)
    write(g, "nzval", v.nzval)
end

function h5_read_sparse(::Type{X}, g::HDF5.Group) where {X<:SparseVector}
    n = read(g["n"])
    nzind = h5_mmap_nonempty(g["nzind"])
    nzval = h5_mmap_nonempty(g["nzval"])
    SparseVector{eltype(nzval),eltype(nzind)}(n, nzind, nzval)
end

function h5_write_sparse(g::HDF5.Group, m::SparseMatrixCSC)
    write(g, "m", m.m)
    write(g, "n", m.n)
    write(g, "colptr", m.colptr)
    write(g, "rowval", m.rowval)
    write(g, "nzval", m.nzval)
end

function h5_read_sparse(::Type{X}, g::HDF5.Group) where {X<:SparseMatrixCSC}
    m = read(g["m"])
    n = read(g["n"])
    colptr = h5_mmap_nonempty(g["colptr"])
    rowval = h5_mmap_nonempty(g["rowval"])
    nzval = h5_mmap_nonempty(g["nzval"])
    SparseMatrixCSC{eltype(nzval),eltype(colptr)}(m, n, colptr, rowval, nzval)
end
