"""
$(README)
"""
module COBREXA

using DocStringExtensions

import AbstractFBCModels as A
import ConstraintTrees as C
import JuMP as J

include("types.jl")
include("solver.jl")

include("builders/core.jl")
include("builders/genes.jl")
include("builders/objectives.jl")

include("analysis/flux_balance_analysis.jl")
include("analysis/parsimonious_flux_balance_analysis.jl")

end # module COBREXA
