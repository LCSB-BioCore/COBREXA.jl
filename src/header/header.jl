
"""
loadSource(srcDirs, root=@__DIR__)

Load source files
"""
function loadSource(srcDirs, rootDir=@__DIR__, testMode=false)
    returnFlag = []
    for d in srcDirs
        srcDir = joinpath(rootDir, d)
        if testMode
            println(" > srcDir = $srcDir")
        end
        for file in filter(f -> endswith(f, ".jl"), readdir(srcDir))
            if file != "header.jl"
                srcFile = joinpath(rootDir, d, file)
                if !testMode
                    include(srcFile)
                else
                    println(" > srcFile = $srcFile")
                end
                push!(returnFlag, true)
            end
        end
    end
    return returnFlag
end