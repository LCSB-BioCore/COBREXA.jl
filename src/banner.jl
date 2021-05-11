const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

function _print_banner()
    c = Base.text_colors
    n = c[:normal]
    y = c[:bold] * c[:cyan]

    println(
        "
                   $(y)//   $(n) |
        \\\\\\\\\\  // $(y)//    $(n) | $(c[:bold])COBREXA.jl $(c[:normal])
         \\\\ \\\\// $(y)//     $(n) |
          \\\\ \\/ $(y)//      $(n) | $(c[:bold])CO$(c[:normal])nstraint-$(c[:bold])B$(c[:normal])ased $(c[:bold])R$(c[:normal])econstruction
           \\\\  $(y)//       $(n) | and $(c[:bold])EX$(c[:normal])ascale $(c[:bold])A$(c[:normal])nalysis in Julia
           //  $(y)\\\\       $(n) |
          // $(y)/\\ \\\\      $(n) | Documentation: https://lcsb-biocore.github.io/COBREXA.jl
         // $(y)//\\\\ \\\\      $(n)|
        // $(y)//  \\\\\\\\\\     $(n)| v$(COBREXA_VERSION) - $(Dates.year(Dates.today()))
       //                $(n)|
        ",
    )
end
