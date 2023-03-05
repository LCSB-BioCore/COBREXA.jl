
"""
    module Types

Module with all types, mainly the model types and various typed model contents.

# Exports
$(EXPORTS)
"""
module Types
using ..ModuleTools
@dse
using ..Internal

using HDF5
using JSON
using LinearAlgebra
using MAT
using OrderedCollections
using SBML
using Serialization
using SparseArrays

@inc_dir types abstract
@export_locals
end # module Types

"""
    module Accessors

Functions that gather data from various model types in a standardized form.
Overload these if you want COBREXA to work with your own module types.

# Exports
$(EXPORTS)
"""
module Accessors
using ..ModuleTools
@dse
using ..Types
using ..Internal.Macros

using SparseArrays

"""
    module Internal

Internal helpers for accessors.

# Exports
$(EXPORTS)
"""
module Internal
    using ..ModuleTools
    @dse
    import ...Types
    import ..Accessors
    using SparseArrays
    @inc_dir types accessors bits
    @export_locals
end

using .Internal

@inc_dir types accessors
@export_locals
end # module Accessors

# the modules depend on each other so we have to inject the stuff like this
@inject Types begin
    using ..Accessors
    using ..Internal.Macros
    using ..Log.Internal: @io_log

    @inc_dir types
    @inc_dir types models

    @export_locals
end

@inject Types.Internal begin
    # TODO where is this declared?
    using ...Types
    using ..Accessors
    using ..Log.Internal: @models_log

    using SBML
    using SparseArrays
    import PikaParser as PP
    using OrderedCollections

    @inc_dir types misc
    @export_locals
end

@inject Types begin
    using .Internal
end
