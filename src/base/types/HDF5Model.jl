
"""
$(TYPEDEF)

A model that is stored in HDF5 format. The model data is never really pulled
into memory, but instead mmap'ed as directly as possible into the Julia
structures.  This makes reading the `HDF5Model`s extremely fast, at the same
time the (uncached) `HDF5Model`s can be sent around efficiently among
distributed nodes just like [`Serialized`](@ref) models, provided the nodes
share a common storage.

All HDF5Models must have the backing disk storage. To create one, use
[`save_h5_model`](@ref) or [`save_model`](@ref) with `.h5` file extension. To
create a temporary model that behaves like a model "in memory", save it to a
temporary file. For related reasons, you can not use `convert` models to
`HDF5Model` format, because the conversion would impliy having the model saved
somewhere.

# Fields
$(TYPEDFIELDS)
"""
mutable struct HDF5Model <: MetabolicModel
    h5::Maybe{HDF5.File}
    filename::String

    HDF5Model(filename::String) = new(nothing, filename)
end
