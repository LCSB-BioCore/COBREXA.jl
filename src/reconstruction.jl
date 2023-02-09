
"""
    module Reconstruction

Reconstruction functionality of COBREXA; mainly functions that construct or
modify the model contents.

# Exports
$(EXPORTS)
"""
module Reconstruction
using ..ModuleTools
@dse

using ..Accessors
using ..Internal: constants
using ..Internal.Macros
using ..Log.Internal
using ..Types

using SparseArrays, OrderedCollections
using MacroTools

@inc_dir reconstruction

"""
    module Pipes

Functions that create model variants, typically for efficient use in
[`screen`](@ref) and similar functions.

# Exports
$(EXPORTS)
"""
module Pipes
    using ..ModuleTools
    @dse
end

@export_locals
end

# this needs to import from Reconstruction
@inject Reconstruction.Pipes begin
    using ..Reconstruction
    using ..Types
    @inc_dir reconstruction pipes

    @export_locals
end

@inject Analysis import ..Reconstruction
