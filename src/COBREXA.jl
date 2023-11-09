"""
$(README)
"""
module COBREXA

using DocStringExtensions

import ConstraintTrees as C
import JuMP as J

include("types.jl")
include("solver.jl")
include("builders/core.jl")
include("builders/objectives.jl")
include("utils/downloads.jl")

end # module COBREXA
