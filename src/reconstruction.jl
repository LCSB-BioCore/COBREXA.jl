
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
using ..Analysis
using ..Internal: constants
using ..Internal.Macros
using ..Log.Internal
using ..Types

using SparseArrays, OrderedCollections
using MacroTools

@inc_dir reconstruction

"""
    module Modifications

Functions that create model variants, typically for efficient use in
[`screen`](@ref) and similar functions.

# Exports
$(EXPORTS)
"""
module Modifications
using ..ModuleTools
@dse
end

@export_locals
end

# this needs to import from Reconstruction
@inject Reconstruction.Modifications begin
    using ...Reconstruction
    @inc_dir reconstruction modifications

    @export_locals
end
