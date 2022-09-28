
# make the module for types and load basic abstract types
module Types
using ..ModuleTools
@dse
using ..Internal
using SparseArrays, OrderedCollections
using HDF5, SBML, JSON, MAT, Serialization #for the storable types
@inc_dir types abstract
@export_locals
end

# the specialized module for accessors
module Accessors
using ..ModuleTools
@dse
using ..Types
using ..Internal.Macros
using SparseArrays

module Internal
using ..ModuleTools
@dse
end

@inc_dir types accessors
@export_locals
end

# the modules depend on each other so we have to inject the stuff like this
@inject Types begin
    using ..Accessors
    using ..Internal.Macros
    using ..Log.Internal: @io_log

    @inc_dir types
    @inc_dir types models
    @inc_dir types wrappers

    @export_locals
end

@inject Types.Internal begin
    using ..Types
    using ..Accessors
    using ..Log.Internal: @models_log

    using SBML
    using SparseArrays

    @inc_dir types misc
    @export_locals
end

@inject Types begin
    using .Internal
end
