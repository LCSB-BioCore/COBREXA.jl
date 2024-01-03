
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
import LinearAlgebra
import SparseArrays

include("types.jl")

# core functionality
include("io.jl")
include("solver.jl")

# conversion of various stuff to constraint trees
include("builders/communities.jl")
include("builders/enzymes.jl")
include("builders/fbc.jl")
include("builders/knockouts.jl")
include("builders/loopless.jl")
include("builders/objectives.jl")
include("builders/thermodynamic.jl")
include("builders/unsigned.jl")
include("builders/variables.jl")

# simplified front-ends for the above
include("frontend/balance.jl")
include("frontend/parsimonious.jl")
include("frontend/mmdf.jl")
include("frontend/moma.jl")
include("frontend/loopless.jl")
include("frontend/enzyme_constrained.jl")

# utilities
include("misc/bounds.jl")
include("misc/settings.jl")
include("misc/utils.jl")

end # module COBREXA
