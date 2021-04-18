
"""
    _formula_to_dict(f::String)::Dict{String,Int}


"""
function _formula_to_dict(f::String)::Dict{String,Int}
    res = Dict{String,Int}()
    element = nothing
    pattern = @r_str "^([A-Z][a-z]*)|([0-9]+)(.*\$)"
    while true
        m = match(pattern, f)

        if isnothing(m)
            break
        elseif !isnothing(m.captures[1])
            if !isnothing(element)
                res[element] = 1
            end
            element = m.captures[1]
        elseif !isnothing(m.captures[2])
            if !isnothing(element)
                res[element] = parse(Int, m.captures[2])
            end
        end
        f = m.captures[3]
    end

    return res
end

"""
    _dict_to_formula(f::String)::Dict{String,Int}


"""
function _dict_to_formula(f::Dict{String,Int})::String
    return join(["$elem$n" for (elem, n) in f])
end
