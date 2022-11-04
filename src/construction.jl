
"""
    module Construction

Model construction functionality of COBREXA; mainly functions that construct or
modify the model contents.

# Exports
$(EXPORTS)
"""
module Construction
using ..ModuleTools
@dse

using ..Accessors
using ..Internal: constants
using ..Internal.Macros
using ..Log.Internal
using ..Types

using SparseArrays, OrderedCollections
using MacroTools

@inc_dir construction

"""
    module With

Functions that create model variants, typically for efficient use in
[`screen`](@ref) and similar functions.

# Exports
$(EXPORTS)
"""
module With
    using ..ModuleTools
    @dse
end

@export_locals
end

# this needs to import from construction
@inject Construction.With begin
    using ..Construction
    @inc_dir construction with_modifications

    @export_locals
end
