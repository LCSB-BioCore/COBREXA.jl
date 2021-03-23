"""
    pfba(model::CobraModel, objective_rxns::Union{Reaction, Array{Reaction, 1}}, optimizer; weights=Float64[], solver_attributes=Dict{Any, Any}(), constraints=Dict{String, Tuple{Float64,Float64}}())

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
    objective_rxns::Union{Reaction,Array{Reaction,1}},
    optimizer;
    weights = Float64[],
    solver_attributes = Dict{Any,Any}(),
    constraints = Dict{String,Tuple{Float64,Float64}}(),
)
    ## FBA ################################################
    cbm, _, _, ubcons, lbcons = build_cbm(model) # get the base constraint based model

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

    # set additional constraints
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, lbcons, ubcons; lb = con[1], ub = con[2])
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
    @objective(cbm, Max, sum(opt_weights[i] * v[i] for i in objective_indices))
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
        0.999999 * λ <= sum(opt_weights[i] * v[i] for i in objective_indices) <= λ
    ) # constrain model - 0.9999 should be close enough?
    @objective(cbm, Min, sum(dot(v, v)))
    optimize!(cbm)

    if termination_status(cbm) != MOI.OPTIMAL &&
       termination_status(cbm) != MOI.LOCALLY_SOLVED # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.99999 * λ <= sum(v[i] for i in objective_indices) <= λ)
        optimize!(cbm)
    end
    if termination_status(cbm) != MOI.OPTIMAL &&
       termination_status(cbm) != MOI.LOCALLY_SOLVED  # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.9999 * λ <= sum(v[i] for i in objective_indices) <= λ)
        optimize!(cbm)
    end
    if termination_status(cbm) != MOI.OPTIMAL &&
       termination_status(cbm) != MOI.LOCALLY_SOLVED  # try to relax bound if failed optimization
        JuMP.delete(cbm, pfbacon)
        @constraint(cbm, 0.999 * λ <= sum(v[i] for i in objective_indices) <= λ)
        optimize!(cbm)
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
