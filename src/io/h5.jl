
"""
$(TYPEDSIGNATURES)

Return a HDF5Model associated with the given file. Does not actually load
anything (for efficiency) -- use [`precache!`](@ref) to start pulling data into
the memory.
"""
function load_h5_model(file_name::String)::HDF5Model
    return HDF5Model(file_name)
end

"""
$(TYPEDSIGNATURES)

Converts and writes a metabolic model to disk in the HDF5 format.

Additionally returns an (uncached) [`HDF5Model`](@ref) that represents the
contents of the saved file. Because all HDF5-based models need to be backed by
disk storage, writing the data to disk (using this function) is the only way to
make new HDF5 models.
"""
function save_h5_model(model::AbstractMetabolicModel, file_name::String)::HDF5Model
    rxns = reactions(model)
    rxnp = sortperm(rxns)
    mets = metabolites(model)
    metp = sortperm(mets)
    h5open(file_name, "w") do f
        write(f, "metabolites", mets[metp])
        write(f, "reactions", rxns[rxnp])
        h5_write_sparse(create_group(f, "balance"), balance(model)[metp])
        h5_write_sparse(create_group(f, "objective"), objective(model)[rxnp])
        h5_write_sparse(create_group(f, "stoichiometry"), stoichiometry(model)[metp, rxnp])
        let (lbs, ubs) = bounds(model)
            write(f, "lower_bounds", lbs[rxnp])
            write(f, "upper_bounds", ubs[rxnp])
        end
    end
    # TODO: genes, grrs, compartments. Perhaps chemistry and others?
    HDF5Model(file_name)
end

"""
$(TYPEDSIGNATURES)

Close (and un-cache) the [`HDF5Model`](@ref) data. This allows the associated
file to be opened for writing again.
"""
function Base.close(model::HDF5Model)
    if !isnothing(model.h5)
        close(model.h5)
        model.h5 = nothing
    end
end
