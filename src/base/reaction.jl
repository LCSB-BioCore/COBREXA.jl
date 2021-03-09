"""
Reaction struct.

# Fields
````
id :: String
name :: String
metabolites :: Dict{Metabolite, Float64}
lb :: Float64
ub :: Float64
grr :: Array{Array{Gene, 1}, 1}
subsystem :: String
notes :: Dict{String, Array{String, 1}}
annotation :: Dict{String, Union{Array{String, 1}, String}}
objective_coefficient :: Float64
````
"""
mutable struct Reaction <: ModelComponent
    id :: String
    name :: String
    metabolites :: Dict{Metabolite, Float64}
    lb :: Float64
    ub :: Float64
    grr :: Array{Array{Gene, 1}, 1}
    subsystem :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}} # everything is a String[] except sbo, which is a String
    objective_coefficient :: Float64
end

"""
Pretty printing of reaction::Reaction.
"""
function Base.show(io::IO, ::MIME"text/plain", r::Reaction)
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

    grr_strings = String[]
    for gr in r.grr
        push!(grr_strings, "("*join([g.id for g in gr], " and ")*")")
    end
    grr_string = join(grr_strings, " or ")
    (isnothing(grr_string) || grr_string == "") && (grr_string = "") 
    println(io, "Genes: ", grr_string)
    println(io, "E.C. number: ", join(get(r.annotation, "ec-code", [""]), " or "))
end

"""
Pretty printing of reactions::Array{Reaction, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", rs::Array{Reaction, 1})
    println(io, "Reaction set of length: ", length(rs))
end

"""
    Reaction()

Empty reaction constructor.

See also: [`Reaction(::String, ::Dict{Metabolite, Float64}, ::String)`](@ref).
"""
function Reaction()
    id = ""
    name = ""
    metabolites = Dict{Metabolite, Float64}()
    lb = -1000.0 
    ub = 1000.0 
    grr = Array{Array{Gene, 1}, 1}()
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

"""
    Reaction(id::String, metabolites::Dict{Metabolite, Float64}, dir="bidir")

Assign the `id`, `metabolites` (and their associated stoichiometries), and the direcionality (`dir`) of a reaction to a `Reaction` struct.
Directionality is specified using `"for"` (forward), `"rev"` (reverse), or any other string for bidirectional reactions. 
All other fields are left unassigned.

See also: [`Reaction()`](@ref).
"""
function Reaction(id::String, metabolites::Dict{Metabolite, Float64}, dir="bidir")
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
    grr = Array{Array{Gene, 1}, 1}()
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

"""
    getindex(rxns::Array{Reaction, 1}, rxn::Reaction)

Get the index of a reaction `rxn` in an array of reactions `rxns`, based in `id`. 
Return -1 if no matches found.

Typically used, `index = rxns[rxn]`.
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

Return the reaction with `rxnid` in `rxns` or else `nothing`.

Typically used: `findfirst(model.rxns, rxnid)`.
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
    check_duplicate_reaction(rxns::Array{Reaction, 1}, rxn::Reaction)

Check if `rxn` already exists in `rxns` but has another id. 
Looks through all the reaction equations of `rxns` and compares metabolite `id`s and their stoichiometric coefficients to those of `rxn`.
If `rxn` has the same reaction equation as another reaction in `rxns`, the return true and the index of the first match. 

See also: [`is_mass_balanced`](@ref)
"""
function check_duplicate_reaction(rxns::Array{Reaction, 1}, crxn::Reaction)
    ceq = Dict{String, Float64}(k.id => v for (k, v) in crxn.metabolites)
    for rxn in rxns
        if rxn.id != crxn.id
            req = Dict{String, Float64}(k.id => v for (k, v) in rxn.metabolites)
            if isempty(setdiff(collect(keys(ceq)), collect(keys(req)))) && isempty(setdiff(collect(keys(req)),collect(keys(ceq)))) # same metabolites
                all([req[k] == v for (k, v) in ceq]) && (return true, rxns[rxn])
            end
        else
            return true, rxns[rxn]
        end
    end
    return false, -1
end

"""
    check_duplicate_annotations(rxns::Array{Reaction, 1}, rxn::Gene)

Determine if a `rxn` is has overlapping annotations in `rxns`.
The annotations checked are: ["bigg.reaction", "biocyc", "ec-code", "kegg.reaction", "metanetx.reaction", "rhea", "sabiork", "seed.reaction"].
Return true and the index of the first hit, otherwise false and -1.
"""
function check_duplicate_annotations(rxns::Array{Reaction, 1}, crxn::Reaction)
    inspect_annotations = ["bigg.reaction", "biocyc", "ec-code", "kegg.reaction", "metanetx.reaction", "rhea", "sabiork", "seed.reaction"]
    for rxn in rxns
        for anno in inspect_annotations
            if any([x in get(crxn.annotation, anno, ["c1"]) for x in get(rxn.annotation, anno, ["c2"])])
                return true, rxns[rxn]
            end
        end
    end
    return false, -1
end

"""
    is_mass_balanced(rxn::Reaction)

Checks if `rxn` is atom balanced. Returns a boolean for whether the reaction is balanced,
and the associated balance of atoms for convenience (useful if not balanced).

See also: [`get_atoms`](@ref), [`check_duplicate_reaction`](@ref)
"""
function is_mass_balanced(rxn::Reaction)
    atom_balances = Dict{String, Float64}() # float here because stoichiometry is not Int64
    for (met, stoich) in rxn.metabolites
        atoms = get_atoms(met)
        isempty(atoms) && continue # ignore blanks
        for (k, v) in atoms
            atom_balances[k] = get(atom_balances, k, 0) + v*stoich
        end
    end

    return all(sum(values(atom_balances)) == 0), atom_balances
end