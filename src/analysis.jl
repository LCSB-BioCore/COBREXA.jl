
module Analysis
using ..ModuleTools
@dse

using ..Types
using ..Accessors
using ..Log.Internal: @models_log

using Distributed, DistributedData
using JuMP

@inc_dir analysis
@inc_dir analysis sampling

module Modifications
using ..ModuleTools
@dse
end

@export_locals
end

# this needs to import from Analysis
@inject Analysis.Modifications begin
    using ...Analysis
    using ...Types
    using ...Accessors

    using JuMP

    @inc_dir analysis modifications

    @export_locals
end
