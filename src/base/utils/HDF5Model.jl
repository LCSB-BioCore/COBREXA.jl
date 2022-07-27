
function h5_write_sparse(g::HDF5.Group, v::SparseVector)
    write(g, "n", v.n)
    write(g, "nzind", v.nzind)
    write(g, "nzval", v.nzval)
end

function h5_read_sparse(::Type{SparseVector}, g::HDF5.Group)
    # TODO mmap
    n = read(g["n"])
    nzind = read(g["nzind"])
    nzval = read(g["nzval"])
    SparseVector{eltype(nzval), eltype(nzind)}(n, nzind, nzval)
end

function h5_write_sparse(g::HDF5.Group, m::SparseMatrixCSC)
    write(g, "m", m.m)
    write(g, "n", m.n)
    write(g, "colptr", m.colptr)
    write(g, "rowval", m.rowval)
    write(g, "nzval", m.nzval)
end

function h5_read_sparse(::Type{SparseMatrixCSC}, g::HDF5.Group)
    # TODO mmap
    m = read(g["m"])
    n = read(g["n"])
    colptr = read(g["colptr"])
    rowval = read(g["rowval"])
    nzval = read(g["nzval"])
    SparseMatrixCSC{eltype(nzval), eltype(colptr)}(m, n, colptr, rowval, nzval)
end
