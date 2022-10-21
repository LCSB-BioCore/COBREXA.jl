
# TODO: It might be much much more useful to actually sort out the utils into
# individual namespaces.
"""
    module Utils

Mostly contains various COBREXA functions and helpers, including internal ones.

# Exports

$(EXPORTS)
"""
module Utils
using ..ModuleTools
@dse

using ..Types
using ..Accessors
using ..Internal.Identifiers
using ..Internal: constants
using ..Internal.Macros
using ..Solver

using HDF5
using JuMP
using LinearAlgebra
using OrderedCollections
using Serialization
using SparseArrays

@inc_dir utils

@export_locals
end

@inject Analysis using ...Utils: objective_bounds
@inject Analysis.Modifications using ...Utils: is_boundary
