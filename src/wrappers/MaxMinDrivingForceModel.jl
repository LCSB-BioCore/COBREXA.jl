"""
$(TYPEDEF)

Return a [`MaxMinDrivingForceModel`](@ref) that can be used to perform max-min
driving force analysis. It is based on the work by Noor, et al., "Pathway
thermodynamics highlights kinetic obstacles in central metabolism.", PLoS
computational biology, 2014.

When [`flux_balance_analysis`](@ref) is called in this type of model, the
max-min driving force algorithm is solved i.e. the objective of the model is to
find the maximum minimum Gibbs free energy of reaction across all reactions in
the model by changing the metabolite concentrations. The variables for max-min
driving force analysis are the actual maximum minimum driving force of the
model, the log metabolite concentrations, and the gibbs free energy reaction
potentials across each reaction. Reaction fluxes are assumed constant.The
optimization problem solved is:
```
max min -ΔᵣG
s.t. ΔᵣG = ΔrG⁰ + R T S' ln(C)
     ΔᵣG ≤ 0
     ln(Cₗ) ≤ ln(C) ≤ ln(Cᵤ)
```
where `ΔrG` are the Gibbs energies dissipated by the reactions, R is the gas
constant, T is the temperature, S is the stoichiometry of the model, and C is
the vector of metabolite concentrations (and their respective lower and upper
bounds).

In case no feasible solution exists, `nothing` is returned.

Since biochemical thermodynamics are assumed, the `proton_ids` and `water_ids`
*need* to be specified so that they can be ignored in the calculations.
Effectively this assumes an aqueous environment at constant pH is used.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct MaxMinDrivingForceModel <: AbstractModelWrapper
    "A dictionary mapping ΔrG⁰ to reactions."
    reaction_standard_gibbs_free_energies::Dict{String,Float64} = Dict{String,Float64}()

    "A cycle-free reference reaction flux solution that is used to set the directions of the reactions. For example, this could be generated this using loopless FBA."
    flux_solution::Dict{String,Float64} = Dict{String,Float64}()

    "Metabolite ids of protons."
    proton_ids::Vector{String} = ["h_c", "h_e"]

    "Metabolite ids of water."
    water_ids::Vector{String} = ["h2o_c", "h2o_e"]

    "A dictionationay mapping metabolite ids to concentrations that are held constant."
    constant_concentrations::Dict{String,Float64} = Dict{String,Float64}()

    "A dictionary mapping metabolite ids to constant concentration ratios in the form `(m1, m2) = r === m1/m2 = r`."
    concentration_ratios::Dict{Tuple{String,String},Float64} =
        Dict{Tuple{String,String},Float64}()

    "Global metabolite concentration lower bound."
    concentration_lb = 1e-9

    "Global metabolite concentration upper bound."
    concentration_ub = 100e-3

    "Thermodynamic temperature."
    T::Float64 = constants.T

    "Real gas constant."
    R::Float64 = constants.R

    "Tolerance use to distinguish flux carrying reactions from zero flux reactions."
    small_flux_tol::Float64 = 1e-6

    "Maximum absolute ΔG bound allowed by a reaction."
    max_dg_bound::Float64 = 1000.0

    "Reaction ids that are ignored internally during thermodynamic calculations. This should include water and proton importers."
    ignore_reaction_ids::Vector{String} = String[]

    "Inner metabolic model calculations are based on."
    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::MaxMinDrivingForceModel) = model.inner

Accessors.variables(model::MaxMinDrivingForceModel) =
    ["mmdf"; "log " .* metabolites(model); "ΔG " .* reaction_ids(model)]

Accessors.n_variables(model::MaxMinDrivingForceModel) =
    1 + n_metabolites(model) + reaction_count(model)

Accessors.metabolite_log_concentration_ids(model::MaxMinDrivingForceModel) =
    "log " .* metabolites(model)
Accessors.metabolite_log_concentration_count(model::MaxMinDrivingForceModel) =
    n_metabolites(model)
Accessors.metabolite_log_concentration_variables(model::MaxMinDrivingForceModel) =
    Dict(mid => Dict(mid => 1.0) for mid in "log " .* metabolites(model))

Accessors.gibbs_free_energy_reaction_ids(model::MaxMinDrivingForceModel) =
    "ΔG " .* reaction_ids(model)
Accessors.gibbs_free_energy_reaction_count(model::MaxMinDrivingForceModel) =
    reaction_count(model)
Accessors.gibbs_free_energy_reaction_variables(model::MaxMinDrivingForceModel) =
    Dict(rid => Dict(rid => 1.0) for rid in "ΔG " .* reaction_ids(model))


Accessors.objective(model::MaxMinDrivingForceModel) =
    [1.0; fill(0.0, n_variables(model) - 1)]

function Accessors.balance(model::MaxMinDrivingForceModel)
    # proton water balance
    num_proton_water = length(model.proton_ids) + length(model.water_ids)
    proton_water_vec = spzeros(num_proton_water)

    # constant concentration balance
    const_conc_vec = log.(collect(values(model.constant_concentrations)))

    # ratio balance
    const_ratio_vec = log.(collect(values(model.concentration_ratios)))

    # give dummy dG0 for reactions that don't have data
    dg0s = [
        get(model.reaction_standard_gibbs_free_energies, rid, 0.0) for
        rid in reaction_ids(model)
    ]

    return [
        proton_water_vec
        const_conc_vec
        const_ratio_vec
        dg0s
    ]
end

function Accessors.stoichiometry(model::MaxMinDrivingForceModel)
    var_ids = Internal.original_variables(model)

    # set proton and water equality constraints
    num_proton_water = length(model.proton_ids) + length(model.water_ids)
    proton_water_mat = spzeros(num_proton_water, n_variables(model))
    idxs = indexin([model.proton_ids; model.water_ids], var_ids)
    for (i, j) in enumerate(idxs)
        isnothing(j) && throw(error("Water or proton ID not found in model."))
        proton_water_mat[i, j] = 1.0
    end

    # constant concentration constraints
    const_conc_mat = spzeros(length(model.constant_concentrations), n_variables(model))
    ids = collect(keys(model.constant_concentrations))
    idxs = indexin(ids, var_ids)
    for (i, j) in enumerate(idxs)
        isnothing(j) &&
            throw(DomainError(ids[j], "Constant metabolite ID not found in model."))
        const_conc_mat[i, j] = 1.0
    end

    # add the relative bounds
    const_ratio_mat = spzeros(length(model.concentration_ratios), n_variables(model))
    for (i, (mid1, mid2)) in enumerate(keys(model.concentration_ratios))
        idxs = indexin([mid1, mid2], var_ids)
        any(isnothing.(idxs)) &&
            throw(DomainError((mid1, mid2), "Metabolite ratio pair not found in model."))
        const_ratio_mat[i, first(idxs)] = 1.0
        const_ratio_mat[i, last(idxs)] = -1.0
    end

    # add ΔG relationships
    dgrs = spdiagm(ones(length(reaction_ids(model))))
    S = stoichiometry(model.inner)
    stoich_mat = -(model.R * model.T) * S'
    dg_mat = [spzeros(reaction_count(model)) stoich_mat dgrs]

    return [
        proton_water_mat
        const_conc_mat
        const_ratio_mat
        dg_mat
    ]
end

function Accessors.bounds(model::MaxMinDrivingForceModel)
    var_ids = Internal.original_variables(model)

    lbs = fill(-model.max_dg_bound, n_variables(model))
    ubs = fill(model.max_dg_bound, n_variables(model))

    # mmdf must be positive for problem to be feasible (it is defined as -ΔG)
    lbs[1] = 0.0
    ubs[1] = 1000.0

    # log concentrations
    lbs[2:(1+n_metabolites(model))] .= log(model.concentration_lb)
    ubs[2:(1+n_metabolites(model))] .= log(model.concentration_ub)

    # need to make special adjustments for the constants
    idxs = indexin([model.proton_ids; model.water_ids], var_ids)
    lbs[idxs] .= -1.0
    ubs[idxs] .= 1.0

    # ΔG for each reaction can be any sign, but that is filled before default

    return (lbs, ubs)
end

function Accessors.coupling(model::MaxMinDrivingForceModel)

    # only constrain reactions that have thermo data
    active_rids = Internal.active_reaction_ids(model)
    idxs = Int.(indexin(active_rids, reaction_ids(model)))

    # thermodynamic sign should correspond to the fluxes
    flux_signs = spzeros(length(idxs), reaction_count(model))
    for (i, j) in enumerate(idxs)
        flux_signs[i, j] = sign(model.flux_solution[reaction_ids(model)[j]])
    end

    neg_dg_mat = [
        spzeros(length(idxs)) spzeros(length(idxs), n_metabolites(model)) flux_signs
    ]

    mmdf_mat = sparse(
        [
            -ones(length(idxs)) spzeros(length(idxs), n_metabolites(model)) -flux_signs
        ],
    )

    return [
        neg_dg_mat
        mmdf_mat
    ]
end

function Accessors.coupling_bounds(model::MaxMinDrivingForceModel)
    n = length(Internal.active_reaction_ids(model))
    neg_dg_lb = fill(-model.max_dg_bound, n)
    neg_dg_ub = fill(0.0, n)

    mmdf_lb = fill(0.0, n)
    mmdf_ub = fill(model.max_dg_bound, n)

    return ([neg_dg_lb; mmdf_lb], [neg_dg_ub; mmdf_ub])
end
