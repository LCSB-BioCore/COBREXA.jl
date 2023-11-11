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
include("solver.jl")

include("builders/core.jl")
include("builders/genes.jl")
include("builders/objectives.jl")
include("builders/enzyme.jl")

include("analysis/modifications/optimizer_settings.jl")
include("analysis/flux_balance_analysis.jl")
include("analysis/parsimonious_flux_balance_analysis.jl")
include("analysis/minimize_metabolic_adjustment_analysis.jl")

end # module COBREXA
