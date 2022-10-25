
"""
    module Analysis

Contains the analysis functions of COBREXA.jl. Typically, these take a
[`AbstractMetabolicModel`](@ref), convert it to the solver represenation and run
various optimization tasks on top of it, such as finding an optimum (e.g. in
[`flux_balance_analysis`](@ref)) or multiple optima (e.g.,
[`flux_variability_analysis`](@ref)).

Functions [`screen`](@ref) and [`screen_optmodel_modifications`](@ref) are
special meta-analyses that apply another specified analysis to multiple
systematically generated versions of the same input model.

# Exports
$(EXPORTS)
"""
module Analysis
using ..ModuleTools
@dse

using ..Accessors
using ..Internal: constants
using ..Log.Internal: @models_log
using ..Solver
using ..Types
using ..Types: _GeckoReactionColumn, _SMomentColumn

using Distributed
using DistributedData
using JuMP
using LinearAlgebra
using Random
using SparseArrays
using StableRNGs

@inc_dir analysis
@inc_dir analysis sampling
@inc_dir analysis reconstruction

"""
    module Modifications

Functions that implement well-defined modifications of the optimization model
before solving. These can be used as `modifications` parameter in e.g.
[`flux_balance_analysis`](@ref) and other analysis functions.

# Exports
$(EXPORTS)
"""
module Modifications
# TODO move this into Solver
using ..ModuleTools
@dse
end

@export_locals
end

# this needs to import from Analysis
@inject Analysis.Modifications begin
    using ...Accessors
    using ...Analysis
    using ...Internal: constants
    using ...Solver
    using ...Types

    using JuMP
    using LinearAlgebra
    using SparseArrays

    @inc_dir analysis modifications

    @export_locals
end
