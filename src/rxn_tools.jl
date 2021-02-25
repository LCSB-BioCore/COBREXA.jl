"""
Reaction struct (mutable)
"""
mutable struct Reaction
    id :: String
    name :: String
    metabolites :: Dict{Metabolite, Float64}
    lb :: Float64
    ub :: Float64
    grr :: String
    subsystem :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}} # SBO is a single string term
    objective_coefficient :: Float64
end

"""
reaction = Reaction()

Empty reaction constructor.
"""
function Reaction()
    id = ""
    name = ""
    metabolites = Dict{Metabolite, Float64}()
    lb = -1000.0 
    ub = 1000.0 
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

"""
reaction = Reaction(metabolite_dict::Dict{Metabolite, Float64}, dir::String)

Assign metabolites and their associated stoichiometries from metabolite_dict to a reaction struct.
Directionality is specified using "for" (forward), "rev" (reverse), or "" for bidirectional reactions. 
All other fields are left unassigned.
"""
function Reaction(metabolites::Dict{Metabolite, Float64}, dir::String)
    id = ""
    name = ""
    if dir == "for"
        lb = 0.0 
        ub = 1000.0     
    elseif dir == "rev"
        lb = -1000.0 
        ub = 0.0     
    else
        lb = -1000.0 
        ub = 1000.0     
    end
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

"""
reaction = Reaction(rxn_dict :: Dict{String, Any}, mets::Array{Metabolite, 1})

Assign a reaction struct using rxn_dict and also check that metabolites in this struct exist in the model.
If not a warning is issued and that metabolite is not added to the reaction.
"""
function Reaction(d :: Dict{String, Any}, mets::Array{Metabolite, 1})
    id = ""
    name = ""
    metabolites = Dict{Metabolite, Float64}()
    lb = -1000.0 
    ub = 1000.0 
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    for (k, v) in d
        if k == "id"
            id = v
        elseif k == "name"
            name = v
        elseif k == "metabolites"
            metabolites = Dict{Metabolite, Float64}() 
            for (kk, vv) in v
                ind = findfirst(x->x.id == kk, mets)
                isnothing(ind) ? (@warn "Metabolite $kk not found in reaction assignment."; continue) : nothing
                metabolites[mets[ind]] = vv
            end 
        elseif k == "lower_bound"
            lb = v
        elseif k == "upper_bound"
            ub = v
        elseif k == "gene_reaction_rule"
            grr = v
        elseif k == "subsystem"
            subsystem = v
        elseif k == "notes"
            notes = Dict{String, Array{String, 1}}(kk=>vv for (kk, vv) in v)
        elseif k == "annotation"
            annotation = Dict{String, Union{Array{String, 1}, String}}()
            for (kk, vv) in v
                if typeof(vv) == String
                    annotation[kk] = vv
                else
                    annotation[kk] = convert(Array{String, 1}, vv)
                end
            end
        elseif k == "objective_coefficient"
            objective_coefficient = v
        else
            @warn "Unrecognized reaction field: $k"
        end
    end
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

"""
index = getindex(model::Model, rxn::Reaction)

Get the index of rxn in model. Return -1 if not found.
"""
function Base.getindex(model::Model, rxn::Reaction)
    return model.rxns[rxn]
end

"""
index = getindex((rxns::Array{Reaction, 1}, rxn::Reaction)

Get the index of a reaction in an array of reactions. Return -1 if not found.
"""
function Base.getindex(rxns::Array{Reaction, 1}, rxn::Reaction)
    for i in eachindex(rxns)
        if rxns[i].id == rxn.id
            return i
        end
    end
    return -1
end

"""
findfirst(rxns::Array{Reaction, 1}, rxnid::String)

Return the reaction with rxnid or else `nothing`. Typically used: findfirst(model.rxns, rxnid)
"""
function Base.findfirst(rxns::Array{Reaction, 1}, rxnid::String)
    for i in eachindex(rxns)
        if rxns[i].id == rxnid
            return rxns[i]
        end
    end
    return nothing
end

"""
Pretty printing of reaction::Reaction.
"""
function Base.show(io::IO, r::Reaction)
    if r.ub > 0.0 && r.lb < 0.0
        arrow = " ⟷  "
    elseif r.ub == 0.0 && r.lb < 0.0
        arrow = " ⟵  "
    elseif r.ub > 0.0 && r.lb == 0.0
        arrow = " ⟶  "
    else
        arrow = " →∣←  " # blocked reaction
    end
    substrates = String[]
    products = String[]
    for (k, v) in r.metabolites
        if v < 0.0
            push!(substrates, string(abs(v))*" "*k.id)
        else
            push!(products, string(abs(v))*" "*k.id)
        end
    end
    isempty(substrates) && (substrates = "∅")
    isempty(products) && (products = "∅")
    
    println(io, "Reaction ID: ", r.id)
    println(io, "Reaction name: ", r.name)
    println(io, "Reaction subsystem: ", r.subsystem)
    if length(substrates) > 5 && length(products) > 5
        sp = substrates[1]*" + ... + "*substrates[end]
        pp = products[1]*" + ... + "*products[end]
        println(io, sp*arrow*pp)
    elseif length(substrates) > 5
        sp = substrates[1]*" + ... + "*substrates[end]
        println(io, sp*arrow*join(products, " + "))
    elseif length(products) > 5
        pp = products[1]*" + ... + "*products[end]
        println(io, join(substrates, " + ")*arrow*pp)    
    else
        println(io, join(substrates, " + ")*arrow*join(products, " + "))
    end
    println(io, "Lower bound: ", r.lb)
    println(io, "Upper bound: ", r.ub)
    println(io, "Genes: ", r.grr)
    println(io, "E.C. number: ", join(get(r.annotation, "ec-code", ["N/A"]), " or "))
end

"""
Pretty printing of reactions::Array{Reaction, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", rs::Array{Reaction, 1})
    println(io, "Reaction set of length: ", length(rs))
end

"""
isbalanced, atom_balances = is_mass_balanced(rxn::Reaction)

Checks if rxn is atom balanced. Returns bool and the associated balance for convenience.
"""
function is_mass_balanced(rxn::Reaction)
    atom_balances = Dict{String, Float64}()
    for (met, stoich) in rxn.metabolites
        atoms = get_atoms(met)
        isempty(atoms) && continue # ignore blanks
        for (k, v) in atoms
            atom_balances[k] = get(atom_balances, k, 0) + v*stoich
        end
    end

    return all(sum(values(atom_balances)) == 0), atom_balances
end