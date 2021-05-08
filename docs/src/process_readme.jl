function process_readme()
    quick_start_section = ""
    add_lines = false
    open(joinpath(@__DIR__, "..", "..", "README.md")) do io
        for ln in eachline(io)
            if !add_lines && startswith(ln, "<!--start-->") # switch on
                add_lines = true
                continue # skip this line
            end
            startswith(ln, "<!--end-->") && break
            if add_lines
                quick_start_section *= "# " * ln * "\n"
            end
        end
    end
    return quick_start_section
end

function from_readme(content)
    content = replace(content, "INSERT_QUICK_START" => process_readme())
    return content
end
