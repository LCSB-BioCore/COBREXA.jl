"""
    module Internal

Internal COBREXA.jl functionality, mostly constants and macros.

# Exports
$(EXPORTS)
"""
module Internal
using ..ModuleTools
@dse

@inc macros
@inc_dir misc ontology
@inc_dir misc

@export_locals
end
