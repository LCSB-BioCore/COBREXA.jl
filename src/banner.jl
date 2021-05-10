
const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

function _print_banner()
    c = Base.text_colors
    n = c[:normal]
    y = c[:bold] * c[:cyan]

    println("       
                 $(y)//   $(n) |
      \\\\\\\\\\  // $(y)//    $(n) |
       \\\\ \\\\// $(y)//     $(n) |    
        \\\\ \\/ $(y)//      $(n) | COnstraint-Based Reconstruction
         \\\\  $(y)//       $(n) | and EXascale Analysis in Julia
         //  $(y)\\\\       $(n) | Version: v$(COBREXA_VERSION)
        // $(y)/\\ \\\\      $(n) | 
       // $(y)//\\\\ \\\\      $(n)|
      // $(y)//  \\\\\\\\\\     $(n)|
     //                $(n)|
            ")
end
