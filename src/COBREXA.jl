
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
import SparseArrays: sparse, findnz
import LinearAlgebra: nullspace

include("types.jl")
include("io.jl")
include("solver.jl")

# these functions build or extend constrainttrees of metabolic models
include("builders/core.jl")
include("builders/genes.jl")
include("builders/objectives.jl")
include("builders/enzymes.jl")
include("builders/thermodynamic.jl")
include("builders/loopless.jl")
include("builders/communities.jl")

# these are the one shot analysis functions
include("frontend/flux_balance_analysis.jl")
include("frontend/parsimonious_flux_balance.jl")
include("frontend/minimization_of_metabolic_adjustment_analysis.jl")
include("frontend/enzyme_constrained_flux_balance_analysis.jl")
include("frontend/loopless_flux_balance_analysis.jl")
include("frontend/max_min_driving_force_analysis.jl")

include("misc/modifications.jl")
include("misc/bounds.jl")
include("misc/utils.jl")

end # module COBREXA
