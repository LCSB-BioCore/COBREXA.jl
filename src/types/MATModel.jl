"""
    struct MATModel

Model struct used when importing models saved in `.mat` format.
All fields in model are imported (i.e. no data loss occurs).
However, not all the fields are used by analysis functions.
"""
struct MATModel <: MetabolicModel
    id::String
    m::Dict{String, Any}
end
