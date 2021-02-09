"""
A struct containing global parameters
"""
mutable struct CobraToolsOptions
    verbose :: Bool
    name_space :: String # kegg, bigg, metanetx, metacyc
end
cto = CobraToolsOptions(true, "bigg")

function Base.show(io::IO, cto::CobraToolsOptions)
    vb = cto.verbose ? "loud" : "quiet"
    println(io, "Output level is ", vb)
    println(io, "Name space used is ", cto.name_space)
end

"""
setverbose(verbose::Bool)

Reduce verbosity (@info and @warn are suppressed) of some functions.
This is mostly useful for unit testing to suppress output (lots of model reading and writing that occurs therein).
"""
function setverbose(verbose)
    if verbose
        cto.verbose = true
    else
        cto.verbose = false
    end
end