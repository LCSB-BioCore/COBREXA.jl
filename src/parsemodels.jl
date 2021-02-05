"""
model = readmodel(file_location)

Reads a model file at file_location and returns a constraint based model.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).
"""
function readmodel(file_location)
    if endswith(file_location, ".json")
        @info "Reading a JSON formatted model..."
    elseif endswith(file_location, ".xml")
        @info "Reading an SBML formatted model..."
    elseif endswith(file_location, ".mat")
        @info "Reading a Matlab formatted model..."
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
        return Model()
    end
end

"""
savemodel(model, file_location)

Save model at location file_location. Infers format from file_location extension.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).
"""
function savemodel(model :: Model, file_location :: String)
    if endswith(file_location, ".json")
        @info "Saving a JSON formatted model..."
    elseif endswith(file_location, ".xml")
        @info "Saving an SBML formatted model..."
    elseif endswith(file_location, ".mat")
        @info "Saving a Matlab formatted model..."
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
    end
end
