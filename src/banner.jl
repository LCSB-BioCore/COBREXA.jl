
const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

function _printBanner()
    c = Base.text_colors
    n = c[:normal] # text
    b = c[:bold] * c[:blue]
    r = c[:bold] * c[:red]
    g = c[:bold] * c[:green]
    m = c[:bold] * c[:magenta]
    println(
        "
  ____ ___  ____  ____  _____$(g)__$(n)  $(r)__$(n)    _     |
 / ___/ _ \\| __ )|  _ \\| ____$(g)\\ \\$(n)$(r)/ /$(n)   / \\    | Constraint-Based Reconstruction
| |  | | | |  _ \\| |_) |  _|  $(g)\\$(n)  $(r)/$(n)   / _ \\   | and EXascale Analysis in Julia
| |__| |_| | |_) |  _ <| |___ $(m)/$(n)  $(b)\\$(n)  / ___ \\  |
 \\____\\___/|____/|_| \\_\\_____$(m)/_/$(n)$(b)\\_\\$(n)/_/   \\_\\ | Version: v$(COBREXA_VERSION)
                                             |")
end
