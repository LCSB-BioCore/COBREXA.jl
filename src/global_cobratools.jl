# Global options for package
mutable struct CobraToolsOptions
    verbose :: Bool
end
cto = CobraToolsOptions(true)

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