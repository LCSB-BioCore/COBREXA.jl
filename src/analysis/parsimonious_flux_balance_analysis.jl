"""
    pfba(model::MetabolicModel, optimizer; modifications, qp_solver, qp_solver_attributes)

Run parsimonious flux balance analysis (pFBA) on the `model`.
Note, the `optimizer` must be set to perform the analysis, any JuMP solver will work.
Optionally, specify problem `modifications` as in [`flux_balance_analysis`](@ref).
Also, `qp_solver` can be set to be different from `optimizer`, where the latter is then the LP optimizer only.
Also note that `qp_solver_attributes` is meant to set the attributes for the `qp_solver`.
Any solver attributes changed in `modifications` will only affect he LP solver.
This function automatically relaxes the constraint that the FBA solution objective matches the pFBA solution. 
This is iteratively relaxed like 1.0, 0.999999, 0.9999, etc. of the bound until 0.99 when the function fails and returns nothing.
Return a solved pFBA model.

# Example
```
optimizer = Gurobi.Optimizer
atts = Dict("OutputFlag" => 0)
model = load_model(StandardModel, "iJO1366.json")
biomass = findfirst(model.reactions, "BIOMASS_Ec_iJO1366_WT_53p95M")
sol = pfba(model, biomass, optimizer; solver_attributes=atts)
```
"""
function parsimonious_flux_balance_analysis(
    model::MetabolicModel,
    optimizer;
    modifications = [(model, opt_model) -> nothing],
    qp_solver = (model, opt_model) -> nothing,
    qp_solver_attributes = [(model, opt_model) -> nothing],
)
    # Run FBA
    opt_model = flux_balance_analysis(model, optimizer; modifications = modifications)
    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing # FBA failed

    # Run pFBA
    λ = objective_value(opt_model)
    old_objective = COBREXA.JuMP.objective_function(opt_model)

    qp_solver(model, opt_model) # change the solver if specified, otherwise default argument does nothing
    if typeof(qp_solver_attributes) <: AbstractVector # many modifications
        for mod in qp_solver_attributes
            mod(model, opt_model)
        end
    else # single modification
        qp_solver_attributes(model, opt_model)
    end

    v = opt_model[:x] # fluxes
    @objective(opt_model, Min, sum(dot(v, v)))

    for lbconval in [1.0, 0.999999, 0.99999, 0.9999, 0.999, 0.99] # sequentially relax bound for stability
        λmin = min(lbconval * λ, λ * 1.0 / lbconval)
        λmax = max(lbconval * λ, λ * 1.0 / lbconval)
        @constraint(
            opt_model,
            pfbacon,
            λmin <= old_objective <= λmax # in case of negative constraints
        )
        optimize!(opt_model)
        COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] &&
            break
        COBREXA.JuMP.delete(opt_model, pfbacon)
    end

    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing # pFBA failed

    return opt_model
end

"""
    parsimonious_flux_balance_analysis_dict(model::MetabolicModel, optimizer; modifications, qp_solver, qp_solver_attributes)

Perform parsimonious flux balance analysis on `model` using `optimizer`. 
Returns a vector of fluxes in the same order as the reactions in `model`. 
Calls [`parsimonious_flux_balance_analysis`](@ref) internally.
"""
function parsimonious_flux_balance_analysis_vec(
    model::MetabolicModel,
    optimizer;
    modifications = [(model, opt_model) -> nothing],
    qp_solver = (model, opt_model) -> nothing,
    qp_solver_attributes = [(model, opt_model) -> nothing],
)
    opt_model = parsimonious_flux_balance_analysis(
        model,
        optimizer;
        modifications = modifications,
        qp_solver = qp_solver,
        qp_solver_attributes = qp_solver_attributes,
    )
    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing

    return value.(opt_model[:x])
end

"""
    parsimonious_flux_balance_analysis_dict(model::MetabolicModel, optimizer; modifications, qp_solver, qp_solver_attributes)

Perform parsimonious flux balance analysis on `model` using `optimizer`. 
Returns a dictionary mapping reaction `id`s to fluxes. 
Calls [`parsimonious_flux_balance_analysis`](@ref) internally.
"""
function parsimonious_flux_balance_analysis_dict(
    model::MetabolicModel,
    optimizer;
    modifications = [(model, opt_model) -> nothing],
    qp_solver = (model, opt_model) -> nothing,
    qp_solver_attributes = [(model, opt_model) -> nothing],
)
    opt_model = parsimonious_flux_balance_analysis(
        model,
        optimizer;
        modifications = modifications,
        qp_solver = qp_solver,
        qp_solver_attributes = qp_solver_attributes,
    )
    COBREXA.JuMP.termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing

    return Dict(zip(reactions(model), value.(opt_model[:x])))
end
