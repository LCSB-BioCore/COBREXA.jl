"""
    module COBREXA

COnstraint Based Reconstruction and EXascale Analysis. COBREXA provides
functions for construction, modification, simulation and analysis of
constraint-based metabolic models that follows the COBRA methodology.

COBREXA is built as a front-end for the combination of `AbstractFBCModels.jl`
(provides the model I/O), `ConstraintTrees.jl` (provides the constraint system
organization), `Distributed.jl` (provides HPC execution capability), and
`JuMP.jl` (provides the solvers).

See the online documentation for a complete description of functionality aided
by copy-pastable examples.

To start quickly, load your favorite JuMP-compatible solver, use
[`load_model`](@ref) to read a metabolic model from the disk, and solve it with
[`flux_balance`](@ref).
"""
module COBREXA

using DocStringExtensions

import AbstractFBCModels as A
import ConstraintTrees as C
import JuMP as J
import SparseArrays: sparse

include("types.jl")
include("io.jl")
include("solver.jl")

include("builders/core.jl")
include("builders/genes.jl")
include("builders/objectives.jl")
include("builders/enzymes.jl")

include("analysis/modifications.jl")
include("analysis/flux_balance.jl")
include("analysis/parsimonious_flux_balance.jl")
include("analysis/minimal_metabolic_adjustment.jl")

include("misc/bounds.jl")

end # module COBREXA
