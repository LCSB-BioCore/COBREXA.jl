module IO
using ..ModuleTools
@dse

using ..Types
using ..Types.Internal: maybemap, unparse_grr
using ..Accessors
using ..Internal: constants
using ..Log.Internal: @io_log

using JSON, MAT, SBML, HDF5

@inc_dir io
@inc_dir io show

module Internal
using ..ModuleTools
@dse
end

@export_locals
end

@inject IO.Internal begin
    using ..Types
    using HDF5
    using SparseArrays

    @inc_dir io misc

    @export_locals
end

@inject IO using .Internal
@inject Types using ..IO.Internal: h5_read_sparse
