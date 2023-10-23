"""
$(README)
"""
module COBREXA

using DocStringExtensions
import AbstractFBCModels as A
import ConstraintTrees as C
import JuMP as J
import SparseArrays: sparse

include("types/maybe.jl")
include("types/isozyme.jl")
include("solver.jl")
include("builders/core.jl")
include("builders/enzyme.jl")

include("utils/downloads.jl")

end # module COBREXA
