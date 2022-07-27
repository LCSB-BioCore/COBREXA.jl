
"""
    load_h5_model(file_name::String)::HDF5Model

Return a HDF5Model associated with the given file. Does not actually load
anything (for efficiency) -- use [`precache!`](@ref) to start pulling data into
the memory.
"""
function load_h5_model(file_name::String)::HDF5Model
    return HDF5Model(file_name)
end

"""
    save_h5_model(model::MetabolicModel, file_name::String)



Additionally returns an (uncached) [`HDF5Model`](@ref) that represents the
contents of the saved file.
"""
function save_h5_model(model::MetabolicModel, file_name::String)::HDF5Model
    h5open(file_name, "w") do f
        write(f, "metabolites", metabolites(model))
        write(f, "reactions", reactions(model))
        h5_write_sparse(create_group(f, "balance"), balance(model))
        h5_write_sparse(create_group(f, "objective"), objective(model))
        h5_write_sparse(create_group(f, "stoichiometry"), stoichiometry(model))
        let (lbs, ubs) = bounds(model)
            write(f, "lower_bounds", lbs)
            write(f, "upper_bounds", ubs)
        end
    end
end

"""
    Base.close(model::HDF5Model)

Close (and un-cache) the [`HDF5Model`](@ref) data. This allows the file to be
opened for writing again.
"""
function Base.close(model::HDF5Model)
    if !isnothing(model.h5)
        close(model.h5)
        model.h5=nothing
    end
end
        
