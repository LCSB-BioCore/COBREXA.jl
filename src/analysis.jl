
module Analysis
using ..ModuleTools
@dse

using ..Accessors
using ..Log.Internal: @models_log
using ..Solver
using ..Types

using Distributed, DistributedData
using JuMP
using StableRNGs, Random

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
    using ...Solver
    using ...Types
    using ...Internal: constants

    using JuMP
    using SparseArrays

    @inc_dir analysis modifications

    @export_locals
end
