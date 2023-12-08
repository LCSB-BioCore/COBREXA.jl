"""
$(README)
"""
module COBREXA

using DocStringExtensions

import AbstractFBCModels as A
import ConstraintTrees as C
import JuMP as J
import SparseArrays: sparse
import LinearAlgebra: dot

include("types.jl")
include("io.jl")
include("solver.jl")

include("builders/core.jl")
include("builders/genes.jl")
include("builders/objectives.jl")

include("analysis/modifications.jl")
include("analysis/flux_balance.jl")
include("analysis/parsimonious_flux_balance.jl")
include("analysis/minimize_metabolic_adjustment.jl")

end # module COBREXA
