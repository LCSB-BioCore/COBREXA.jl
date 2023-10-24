"""
$(README)
"""
module COBREXA

using DocStringExtensions

import ConstraintTrees as C

include("types.jl")
include("solver.jl")
include("builders/core.jl") #TODO more stuff
include("utils/downloads.jl")

end # module COBREXA
