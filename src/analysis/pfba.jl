"""
    pfba(model::CobraModel, optimizer; objective_func::Union{Reaction, Array{Reaction, 1}}=Reaction[], weights=Float64[], solver_attributes=Dict{Any, Any}(), constraints=Dict{String, Tuple{Float64,Float64}}(), sense=MOI.MAX_SENSE)

Run parsimonious flux balance analysis (pFBA) on the `model` with `objective_rxn(s)` and optionally specifying their `weights` (empty `weights` mean equal weighting per reaction) for the initial FBA problem.
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work.
Optionally also specify any additional flux constraints with `constraints`, a dictionary mapping reaction `id`s to tuples of (lb, ub) flux constraints.
When `optimizer` is an array of optimizers, e.g. `[opt1, opt2]`, then `opt1` is used to solve the FBA problem, and `opt2` is used to solve the QP problem.
This strategy is useful when the QP solver is not good at solving the LP problem.
The `solver_attributes` can also be specified in the form of a dictionary where each (key, value) pair will be passed to `set_optimizer_attribute(cbmodel, k, v)`.
If more than one solver is specified in `optimizer`, then `solver_attributes` must be a dictionary of dictionaries with keys "opt1" and "opt2", e.g. Dict("opt1" => Dict{Any, Any}(),"opt2" => Dict{Any, Any}()).
This function builds the optimization problem from the model, and hence uses the constraints implied by the model object.
Returns a dictionary of reaction `id`s mapped to fluxes if solved successfully, otherwise an empty dictionary.

# Example
```
optimizer = Gurobi.Optimizer
atts = Dict("OutputFlag" => 0)
model = CobraTools.read_model("iJO1366.json")
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")
sol = pfba(model, biomass, optimizer; solver_attributes=atts)
```
"""
function pfba(
    model::CobraModel,
    optimizer;
    modifications = [(model, opt_model) -> nothing]
)
    ## FBA ################################################

    if typeof(optimizer) <: AbstractArray # choose optimizer
        cbm = make_optimization_model(model, optimizer[1], sense = sense)
        v = cbm[:x]
        if !isempty(solver_attributes["opt1"]) # set other attributes
            for (k, v) in solver_attributes["opt1"]
                set_optimizer_attribute(cbm, k, v)
            end
        end
    else # singe optimizer
        cbm = make_optimization_model(model, optimizer, sense = sense)
        v = cbm[:x]
        if !isempty(solver_attributes) # set other attributes
            for (k, v) in solver_attributes
                set_optimizer_attribute(cbm, k, v)
            end
        end
    end

    # set additional constraints
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, cbm; lb = con[1], ub = con[2])
    end

    # check if default objective should be used
    if typeof(objective_func) == Reaction || !isempty(objective_func)
        # check if an array of objective indices are fed in
        if typeof(objective_func) == Reaction
            objective_indices = [model[objective_func]]
        else
            objective_indices = [model[rxn] for rxn in objective_func]
        end

        if isempty(weights)
            weights = ones(length(objective_indices))
        end
        opt_weights = zeros(length(model.reactions))

        wcounter = 1
        for i in eachindex(model.reactions)
            if i in objective_indices
                # model.reactions[i].objective_coefficient = weights[wcounter]
                opt_weights[i] = weights[wcounter]
                wcounter += 1
                # else
                # model.reactions[i].objective_coefficient = 0.0
            end
        end

        @objective(cbm, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    else
        # objective_indices = findnz(objective(model))
        # opt_weights = ones(length(objective_indices)) # assume equal weighting, assume sense is max
    end

    optimize!(cbm)

    fba_status = (
        termination_status(cbm) == MOI.OPTIMAL ||
        termination_status(cbm) == MOI.LOCALLY_SOLVED
    )

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

    @constraint(
        cbm,
        pfbacon,
        λ <= sum(opt_weights[i] * v[i] for i in objective_indices) <= λ
    )
    @objective(cbm, Min, sum(dot(v, v)))

    optimize!(cbm)

    for lbconval in [0.999999, 0.99999, 0.9999, 0.999, 0.99] # relax bound for stability
        if termination_status(cbm) == MOI.OPTIMAL ||
           termination_status(cbm) == MOI.LOCALLY_SOLVED # try to relax bound if failed optimization
            break
        else
            COBREXA.JuMP.delete(cbm, pfbacon)
            @constraint(cbm, lbconval * λ <= sum(v[i] for i in objective_indices) <= λ)
            optimize!(cbm)
        end
    end

    pfba_status = (
        termination_status(cbm) == MOI.OPTIMAL ||
        termination_status(cbm) == MOI.LOCALLY_SOLVED
    )

    if fba_status && pfba_status
        return map_fluxes(v, model)
    else
        @warn "Optimization issues occurred."
        return Dict{String,Float64}() # something went wrong
    end
end
