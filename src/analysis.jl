
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

module Modifications
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
