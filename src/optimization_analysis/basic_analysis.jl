"""
    get_core_model(model::CobraTools.Model)

Return stoichiometrix matrix (S), mass balance right hand side (b), upper (ubs), and lower bounds (lbs) of constraint based `model`.
That is, S*v=b with lbs ≤ v ≤ ubs where v is the flux vector. 
This is useful if you want to construct your own optimization problem.

Returns: S, b, ubs, lbs
"""
function get_core_model(model::CobraTools.Model)
    ubs = [rxn.ub for rxn in model.reactions]
    lbs = [rxn.lb for rxn in model.reactions]
    
    b = spzeros(length(model.metabolites))
    S = spzeros(length(model.metabolites), length(model.reactions))

    metids = [met.id for met in model.metabolites] # need indices for S matrix construction
    for (i, rxn) in enumerate(model.reactions) # column
        for (met, coeff) in rxn.metabolites
            j = findfirst(x -> x == met.id, metids) # row
            isnothing(j) ? (@error "S matrix construction error: $(met.id) not defined."; continue) : nothing
            S[j, i] = coeff
        end
    end
    return S, b, ubs, lbs
end

"""
    build_cbm(model::CobraTools.Model)

Initialize a constraint based `model` using `JuMP`. 
Creates a model that satisfies the mass balance and flux constraints but no objective or optimizer is set. 
Returns the `JuMP` model.
This is useful if you want to write your own optimization problem.

Returns: `cbmodel`, `v`, `mb`, `ubs`, and `lbs`, where `cbmodel` is the JuMP model, `v` are the fluxes, `mb` is S*v == 0, and lbs <= v <= ubs.
"""
function build_cbm(model::CobraTools.Model)
    S, b, ubs, lbs = get_core_model(model) # Construct S, b, lbs, ubs from model
    cbmodel = JuMP.Model()
    nvars = size(S, 2) # number of variables in model
    v = @variable(cbmodel, v[1:nvars]) # flux variables
    mb = @constraint(cbmodel, mb, S*v .== b) # mass balance
    lbs = @constraint(cbmodel, lbs, lbs .<= v) # lower bounds
    ubs = @constraint(cbmodel, ubs, v .<= ubs) # upper bounds
    return cbmodel, v, mb, ubs, lbs
end

@doc raw"""
    fba(model::CobraTools.Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}, optimizer; weights=Float64[], solver_attributes=Dict{Any, Any}())

Run flux balance analysis (FBA) on the `model` with `objective_rxn(s)` and optionally specifying their `weights` (empty `weights` mean equal weighting per reaction).
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work. 
The `solver_attributes` can also be specified in the form of a dictionary where each (key, value) pair will be passed to `set_optimizer_attribute(cbmodel, key, value)`.
This function builds the optimization problem from the model, and hence uses the constraints implied by the model object.
Returns a dictionary of reaction `id`s mapped to fluxes if solved successfully, otherwise an empty dictionary.

# Example
optimizer = Gurobi.Optimizer

atts = Dict("OutputFlag" => 0)

model = CobraTools.read_model("iJO1366.json")

biomass = findfirst(model.reactions, "BIOMASS\_Ec\_iJO1366\_WT\_53p95M")

sol = fba(model, biomass, optimizer; solver_attributes=atts)
"""
function fba(model::CobraTools.Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}, optimizer; weights=Float64[], solver_attributes=Dict{Any, Any}())
    cbm, _, _, _, _ = build_cbm(model) # get the base constraint based model
    set_optimizer(cbm, optimizer) # choose optimizer
    if !isempty(solver_attributes) # set other attributes
        for (k, v) in solver_attributes
            set_optimizer_attribute(cbm, k, v)
       end
    end

    # ensure that an array of objective indices are fed in
    if typeof(objective_rxns) == Reaction
        objective_indices = [model[objective_rxns]]
    else
        objective_indices = [model[rxn] for rxn in objective_rxns]
    end
    
    if isempty(weights)
        weights = ones(length(objective_indices))
    end
    opt_weights = zeros(length(model.reactions))
    
    # update the objective function tracker
    wcounter = 1
    for i in eachindex(model.reactions)
        if i in objective_indices
            model.reactions[i].objective_coefficient = weights[wcounter]
            opt_weights[i] = weights[wcounter]
            wcounter += 1
        else
            model.reactions[i].objective_coefficient = 0.0
        end
    end
    
    v = all_variables(cbm)
    @objective(cbm, Max, sum(opt_weights[i]*v[i] for i in objective_indices))
    optimize!(cbm)    

    status = (termination_status(cbm) == MOI.OPTIMAL || termination_status(cbm) == MOI.LOCALLY_SOLVED)
    
    if status
        return map_fluxes(v, model)      
    else
        @warn "Optimization issues occurred."
        return Dict{String, Float64}()
    end
end

@doc raw"""
    pfba(model::CobraTools.Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}, optimizer; weights=Float64[], solver_attributes=Dict{Any, Any}())

Run parsimonious flux balance analysis (pFBA) on the `model` with `objective_rxn(s)` and optionally specifying their `weights` (empty `weights` mean equal weighting per reaction) for the initial FBA problem.
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work.
When `optimizer` is an array of optimizers, e.g. `[opt1, opt2]`, then `opt1` is used to solve the FBA problem, and `opt2` is used to solve the QP problem.
This strategy is useful when the QP solver is not good at solving the LP problem.
The `solver_attributes` can also be specified in the form of a dictionary where each (key, value) pair will be passed to `set_optimizer_attribute(cbmodel, k, v)`.
If more than one solver is specified in `optimizer`, then `solver_attributes` must be a dictionary of dictionaries with keys "opt1" and "opt2", e.g. Dict("opt1" => Dict{Any, Any}(),"opt2" => Dict{Any, Any}()).
This function builds the optimization problem from the model, and hence uses the constraints implied by the model object.
Returns a dictionary of reaction `id`s mapped to fluxes if solved successfully, otherwise an empty dictionary.

# Example
optimizer = Gurobi.Optimizer

atts = Dict("OutputFlag" => 0)

model = CobraTools.read_model("iJO1366.json")

biomass = findfirst(model.reactions, "BIOMASS\_Ec\_iJO1366\_WT\_53p95M")

sol = pfba(model, biomass, optimizer; solver_attributes=atts)
"""
function pfba(model::CobraTools.Model, objective_rxns::Union{Reaction, Array{Reaction, 1}}, optimizer; weights=Float64[], solver_attributes=Dict{Any, Any}())
    ## FBA ################################################
    cbm, _, _, _, _ = build_cbm(model) # get the base constraint based model
    
    if typeof(optimizer) <: AbstractArray # choose optimizer
        set_optimizer(cbm, optimizer[1]) 
        if !isempty(solver_attributes["opt1"]) # set other attributes
            for (k, v) in solver_attributes["opt1"]
                set_optimizer_attribute(cbm, k, v)
           end
        end
    else
        set_optimizer(cbm, optimizer) # choose optimizer
        if !isempty(solver_attributes) # set other attributes
            for (k, v) in solver_attributes
                set_optimizer_attribute(cbm, k, v)
           end
        end
    end

    # ensure that an array of objective indices are fed in
    if typeof(objective_rxns) == Reaction
        objective_indices = [model[objective_rxns]]
    else
        objective_indices = [model[rxn] for rxn in objective_rxns]
    end
    
    if isempty(weights)
        weights = ones(length(objective_indices))
    end
    opt_weights = zeros(length(model.reactions))
    
    # update the objective function tracker
    wcounter = 1
    for i in eachindex(model.reactions)
        if i in objective_indices
            model.reactions[i].objective_coefficient = weights[wcounter]
            opt_weights[i] = weights[wcounter]
            wcounter += 1
        else
            model.reactions[i].objective_coefficient = 0.0
        end
    end
    
    v = all_variables(cbm)
    @objective(cbm, Max, sum(opt_weights[i]*v[i] for i in objective_indices))
    optimize!(cbm)    

    fba_status = (termination_status(cbm) == MOI.OPTIMAL || termination_status(cbm) == MOI.LOCALLY_SOLVED)

    ## pFBA ###############################################
    λ = objective_value(cbm)

    if typeof(optimizer) <: AbstractArray # choose optimizer
        set_optimizer(cbm, optimizer[2]) 
        if !isempty(solver_attributes["opt2"]) # set other attributes
            for (k, v) in solver_attributes["opt2"]
                set_optimizer_attribute(cbm, k, v)
           end
        end
    end

    @constraint(cbm, pfbacon, 0.999999*λ <= sum(opt_weights[i]*v[i] for i in objective_indices) <= λ) # constrain model - 0.9999 should be close enough?
    @objective(cbm, Min, sum(dot(v, v)))
    optimize!(cbm)

    if termination_status(cbm) != MOI.OPTIMAL && termination_status(cbm) != MOI.LOCALLY_SOLVED # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.99999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end
    if termination_status(cbm) != MOI.OPTIMAL && termination_status(cbm) != MOI.LOCALLY_SOLVED  # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.9999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end
    if termination_status(cbm) != MOI.OPTIMAL && termination_status(cbm) != MOI.LOCALLY_SOLVED  # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.999*λ <= sum(v[i] for i in objective_indices) <= λ) 
        optimize!(cbm)    
    end

    pfba_status = (termination_status(cbm) == MOI.OPTIMAL || termination_status(cbm) == MOI.LOCALLY_SOLVED)
    
    if fba_status && pfba_status
        return map_fluxes(v, model)      
    else
        @warn "Optimization issues occurred."
        return Dict{String, Float64}() # something went wrong
    end
end

"""
    map_fluxes(v, model::CobraTools.Model)

Map fluxes from an optimization problem (`v`) to rxns in a model. 
`v` can be a JuMP object (fluxes) or an array of Float64 fluxes.
Assumes they are in order of `model.reactions`, which they should be since the optimization problem is constructed from the model.
"""
function map_fluxes(v::Array{VariableRef,1}, model::CobraTools.Model)
    rxndict = Dict{String, Float64}()
    for i in eachindex(model.reactions)
        rxndict[model.reactions[i].id] = value(v[i])
    end
    return rxndict
end

function map_fluxes(v::Array{Float64,1}, model::CobraTools.Model)
    rxndict = Dict{String, Float64}()
    for i in eachindex(model.reactions)
        rxndict[model.reactions[i].id] = v[i]
    end
    return rxndict
end

"""
    setbound(index, ubconstaintref, lbconstaintref; ub=1000, lb=-1000)

Helper function to set the bounds of variables.
The JuMP set_normalized_rhs function is a little confusing, so this function simplifies setting
constraints. Just supply the constraint `index` and the bound objects (`ubconstaintref`, `lbconstaintref`) and they will be set to `ub` and `lb`.

See also: [`set_normalized_rhs`](@ref)
"""
function set_bound(vind, ubs, lbs; ub=1000, lb=-1000)
    if lb <= 0 
        set_normalized_rhs(lbs[vind], abs(lb))
    else
        set_normalized_rhs(lbs[vind], -abs(lb))
    end
    set_normalized_rhs(ubs[vind], ub)
end

"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::CobraTools.Model)

Return a dictionary mapping the flux of atoms across the boundary of the model given `flux_dict` of reactions in `model`. 
Here `flux_dict` is a mapping of reaction `id`s to fluxes, e.g. from FBA.
"""
function atom_exchange(flux_dict::Dict{String, Float64}, model::CobraTools.Model)
    atom_flux = Dict{String, Float64}()
    for (rxnid, flux) in flux_dict
        if startswith(rxnid, "EX_") || startswith(rxnid, "DM_") # exchange, demand reaction
            for (met, stoich) in findfirst(model.reactions, rxnid).metabolites
                adict = get_atoms(met)
                for (atom, stoich) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux*stoich
                end
            end
        end
    end
    return atom_flux
end

"""
    get_exchanges(rxndict::Dict{String, Float64}; topN=8, ignorebound=1000.0, verbose=true)

Display the topN producing and consuming exchange fluxes. 
Set topN to a large number to get all the consuming/producing fluxes.
Ignores infinite (problem upper/lower bound) fluxes (set with ignorebound).
When `verbose` is false, the output is not printed out.
Return these reactions in two dictionaries: `consuming`, `producing`
"""
function exchange_reactions(rxndict::Dict{String, Float64}; topN=8, ignorebound=1000.0, verbose=true)
    fluxes = Float64[]
    rxns = String[]
    for (k, v) in rxndict
        if startswith(k, "EX_") && abs(v) < ignorebound
            push!(rxns, k)
            push!(fluxes, v)
        end
    end
    inds_prod = sortperm(fluxes, rev=true)
    inds_cons = sortperm(fluxes)

    consuming = Dict{String, Float64}()
    producing = Dict{String, Float64}()
    verbose && println("Consuming fluxes:")
    for i in 1:min(topN, length(rxndict))
        if rxndict[rxns[inds_cons[i]]] < -eps()
            verbose && println(rxns[inds_cons[i]], " = ", round(rxndict[rxns[inds_cons[i]]], digits=4))
            consuming[rxns[inds_cons[i]]] = rxndict[rxns[inds_cons[i]]]
        else
            continue
        end
    end

    verbose && println("Producing fluxes:")
    for i in 1:min(topN, length(rxndict))
        if rxndict[rxns[inds_prod[i]]] > eps()
            verbose && println(rxns[inds_prod[i]], " = ", round(rxndict[rxns[inds_prod[i]]], digits=4))
            producing[rxns[inds_prod[i]]] = rxndict[rxns[inds_prod[i]]]
        else
            continue
        end
    end
    return consuming, producing
end

"""
    metabolite_fluxes(fluxdict::Dict{String, Float64}, model::CobraTools.Model)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or produce them.
"""
function metabolite_fluxes(fluxdict::Dict{String, Float64}, model::CobraTools.Model)
    S, _, _, _ = get_core_model(model)
    S = Array(S) # full
    met_flux = Dict{String, Float64}()
    rxnids = [rxn.id for rxn in model.reactions]
    metids= [met.id for met in model.metabolites]

    producing = Dict{String, Dict{String, Float64}}()
    consuming = Dict{String, Dict{String, Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = fluxdict[rxnid]*S[row, col]
            # ignore zero flux
            if mf < -eps() # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > eps()
                if haskey(producing, metid)
                    producing[metid][rxnid] = mf
                else
                    producing[metid] = Dict(rxnid => mf)                    
                end                
            end
        end
    end
    return consuming, producing
end
