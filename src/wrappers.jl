
"""
    module Wrappers

All "layered" modifications of the models and their types are grouped in this module.

# Exports
$(EXPORTS)
"""
module Wrappers
using ..ModuleTools
@dse

using ..Types
using ..Accessors
using ..Internal.Macros
using ..Internal: constants

using LinearAlgebra
using SparseArrays

module Internal
    using ..ModuleTools
    @dse
    using ..Wrappers
    using ..Types
    using ..Accessors

    using SparseArrays

    @inc_dir wrappers bits
    @export_locals
end # module Internal

using .Internal

@inc_dir wrappers

@export_locals
end # module Wrappers

@inject Wrappers.Internal begin
    @inc_dir wrappers misc
    @export_locals # again
end
