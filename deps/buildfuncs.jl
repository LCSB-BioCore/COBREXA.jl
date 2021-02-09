"""
mkGibbsDB(equilibrator_file_loc::String)

Read and parse equilibrator ΔG data. Save to Julia file. 
"""
function mkGibbsDB(in_file_loc::String, out_file_loc::String)
    gibbs = Dict{String, Measurement{Float64}}()
    open(in_file_loc) do io
        for ln in eachline(io)
            startswith(ln, "!") && continue
            prts = split(ln, "\t")
            if prts[2] != "nan"
                gibbs[prts[1]] = parse(Float64, prts[2]) ± parse(Float64, prts[3]) 
            end
        end
    end
    JLD.save(out_file_loc, "gibbs", gibbs)
end
