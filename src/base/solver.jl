
"""
    make_optimization_model(
        model::MetabolicModel,
        optimizer;
        sense = MOI.MAX_SENSE,
    )

Convert `MetabolicModel`s to a JuMP model, place objectives and the equality
constraint.
"""
function make_optimization_model(model::MetabolicModel, optimizer; sense = MOI.MAX_SENSE)

    precache!(model)

    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = Model(optimizer)
    @variable(optimization_model, x[i = 1:n])
    @objective(optimization_model, sense, objective(model)' * x)
    @constraint(optimization_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    @constraint(optimization_model, lbs, xl .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= xu) # upper bounds

    C = coupling(model) # empty if no coupling
    cl, cu = coupling_bounds(model)
    isempty(C) || @constraint(optimization_model, c_lbs, cl .<= coupling(model) * x) # coupling lower bounds
    isempty(C) || @constraint(optimization_model, c_ubs, coupling(model) * x .<= cu) # coupling upper bounds

    enzyme_vec, enzyme_mass = enzyme_capacity(model) # nothing if not present
    !isnothing(enzyme_vec) &&
        @constraint(optimization_model, enz_cap, dot(enzyme_vec, x) <= enzyme_mass)

    return optimization_model
end

"""
    is_solved(opt_model)

Return `true` if `opt_model` solved successfully (solution is optimal or locally
optimal).  Return `false` if any other termination status is reached.
Termination status is defined in the documentation of `JuMP`.
"""
is_solved(opt_model) = termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]

"""
    optimize_objective(opt_model)::Maybe{Float64}

Shortcut for running JuMP `optimize!` on a model and returning the objective
value, if solved.
"""
function optimize_objective(opt_model)::Maybe{Float64}
    optimize!(opt_model)
    solved_objective_value(opt_model)
end

"""
    get_optmodel_bounds(opt_model)

Returns vectors of the lower and upper bounds of `opt_model` constraints, where
`opt_model` is a JuMP model constructed by e.g.
[`make_optimization_model`](@ref) or [`flux_balance_analysis`](@ref).
"""
get_optmodel_bounds(opt_model) = (
    [-normalized_rhs(lb) for lb in opt_model[:lbs]],
    [normalized_rhs(ub) for ub in opt_model[:ubs]],
)

"""
    set_optmodel_bound!(vidx, opt_model;
        ub::Maybe{Real} = nothing,
        lb::Maybe{Real} = nothing,
    )

Helper function to set the bounds of a variable in the model. Internally calls
`set_normalized_rhs` from JuMP. If the bounds are set to `nothing`, they will
not be changed.
"""
function set_optmodel_bound!(
    vidx,
    opt_model;
    lb::Maybe{Real} = nothing,
    ub::Maybe{Real} = nothing,
)
    isnothing(lb) || set_normalized_rhs(opt_model[:lbs][vidx], -lb)
    isnothing(ub) || set_normalized_rhs(opt_model[:ubs][vidx], ub)
end


"""
    solved_objective_value(opt_model)::Maybe{Float64}

Returns the current objective value of a model, if solved.

# Example
```
solved_objective_value(flux_balance_analysis(model, ...))
```
"""
solved_objective_value(opt_model)::Maybe{Float64} =
    is_solved(opt_model) ? objective_value(opt_model) : nothing

"""
    flux_vector(opt_model)::Maybe{Vector{Float64}}

Returns a vector of fluxes of the model, if solved.

# Example
```
flux_vector(flux_balance_analysis(model, ...))
```
"""
flux_vector(model::MetabolicModel, opt_model)::Maybe{Vector{Float64}} =
    is_solved(opt_model) ? reaction_flux(model)' * value.(opt_model[:x]) : nothing

"""
    flux_dict(model::MetabolicModel, opt_model)::Maybe{Dict{String, Float64}, Nothing}

Returns the fluxes of the model as a reaction-keyed dictionary, if solved.

# Example
```
flux_dict(model, flux_balance_analysis(model, ...))
```
"""
flux_dict(model::MetabolicModel, opt_model)::Maybe{Dict{String,Float64}} =
    is_solved(opt_model) ?
    Dict(reactions(model) .=> reaction_flux(model)' * value.(opt_model[:x])) : nothing

"""
    flux_dict(model::GeckoModel, opt_model)

Specialization to format solved data for `GeckoModel`s but maps 
the solution back into the namespace of the underlying model (the 
original ids).
"""
flux_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    _map_irrev_to_rev_ids(model.geckodata.reaction_map, value.(opt_model[:x])) : nothing


"""
    flux_dict(model::SMomentModel, opt_model)

Specialization to format solved data for `SMomentModel`s but maps 
the solution back into the namespace of the underlying model (the 
original ids).
"""
flux_dict(model::SMomentModel, opt_model) =
    is_solved(opt_model) ?
    _map_irrev_to_rev_ids(model.smomentdata.reaction_map, value.(opt_model[:x])) : nothing

"""
    _map_irrev_to_rev_ids(reaction_map, protein_ids, solution)

Return dictionaries of reaction ids mapped to fluxes, 
and protein ids mapped to concentrations using `reaction_map` to 
determine the ids of fluxes and `protein_ids` for the gene ids.
The solution in `solution` is used to fill the dictionaries.
"""
function _map_irrev_to_rev_ids(reaction_map, solution; protein_ids = [])
    reaction_flux = Dict{String,Float64}()
    for (k, i) in reaction_map
        contains(k, "§ISO") && continue # §ISO§FOR and §ISO§REV need to be ignored
        rid = split(k, "§")[1]
        v = contains(k, "§FOR") ? solution[i] : -solution[i]
        reaction_flux[rid] = get(reaction_flux, rid, 0) + v
    end

    if isempty(protein_ids)
        return reaction_flux
    else
        n_reactions = length(reaction_map)
        protein_flux = Dict{String,Float64}()
        for (i, pid) in enumerate(protein_ids)
            protein_flux[pid] = solution[n_reactions+i]
        end
        return reaction_flux, protein_flux
    end
end
