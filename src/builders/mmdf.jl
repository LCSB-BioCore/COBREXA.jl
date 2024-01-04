
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
$(TYPEDSIGNATURES)

Build a max min driving force analysis model. See the docstring of
[`max_min_driving_force_analysis`](@ref) for details about the arguments.
"""
function max_min_driving_force_constraints(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    reference_flux = Dict{String,Float64}(),
    proton_ids::Vector{String},
    water_ids::Vector{String},
    concentration_lb = 1e-9, # M
    concentration_ub = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = String[],
)
    # check if reactions that will be used in the model all have thermodynamic data, otherwise throw error.
    # only these reactions will form part of the model.
    rxn_source =
        isempty(reference_flux) ? A.reactions(model) : collect(keys(reference_flux))
    rxns = Set(x for x in rxn_source if !(x in ignore_reaction_ids))
    isempty(setdiff(rxns, keys(reaction_standard_gibbs_free_energies))) || throw(
        DomainError(
            reaction_standard_gibbs_free_energies,
            """
                      Not all reactions have thermodynamic data.
                      Either add the reactions with missing ΔG0s to ignore_reaction_ids,
                      or add the data to reaction_standard_gibbs_free_energies.
                      """,
        ),
    )
    rxns = collect(rxns)

    # import metabolite ids (if used), and reaction stoichiometries from an AbstractFBCModel
    mets = Set(m for rid in rxns for m in keys(A.reaction_stoichiometry(model, rid)))
    mets = collect(mets)
    stoi = A.stoichiometry(model)

    # create thermodynamic variables
    constraints = C.ConstraintTree(
        :max_min_driving_force^C.variable() +
        :log_metabolite_concentrations^C.variables(
            keys = Symbol.(mets),
            bounds = C.Between(log(concentration_lb), log(concentration_ub)),
        ) +
        :delta_G_reactions^C.variables(keys = Symbol.(rxns)),
    )

    #=
    Build gibbs free energy relations
    ΔG_rxn == ΔG0 + R * T * sum ( log_concentration_variable * stoichiometry_value )
    =#
    model_rxns = A.reactions(model)
    model_mets = A.metabolites(model)
    dG0s_met_ids_stoichs = Vector{Tuple{Float64,Vector{String},Vector{Float64}}}()
    for rxn in rxns # prepare to create CT below
        met_idxs, stoich_coeffs =
            SparseArrays.findnz(stoi[:, findfirst(==(rxn), model_rxns)])
        met_ids = model_mets[met_idxs]
        dG0 = reaction_standard_gibbs_free_energies[rxn]
        push!(dG0s_met_ids_stoichs, (dG0, met_ids, stoich_coeffs))
    end

    constraints *=
        :delta_G_reaction_equations^C.ConstraintTree(
            Symbol(rxn) => C.Constraint(
                value = -constraints.delta_G_reactions[Symbol(rxn)].value +
                        dG0 +
                        R *
                        T *
                        sum(
                            constraints.log_metabolite_concentrations[Symbol(
                                met_id,
                            )].value * stoich for
                            (met_id, stoich) in zip(met_ids, stoich_coeffs)
                        ),
                bound = C.EqualTo(0.0),
            ) for (rxn, (dG0, met_ids, stoich_coeffs)) in zip(rxns, dG0s_met_ids_stoichs)
        )

    #=
    Set proton log concentration to zero so that it won't impact any
    calculations (biothermodynamics assumption). Also set water concentrations
    to zero (aqueous condition assumptions). How true is "aqueous conditions"?
    Debatable...
    =#
    for met in [Symbol.(proton_ids); Symbol.(water_ids)]
        if haskey(constraints.log_metabolite_concentrations, met)
            constraints.log_metabolite_concentrations[met] = C.Constraint(
                value = constraints.log_metabolite_concentrations[met].value,
                bound = C.EqualTo(0.0),
            )
        end
    end

    #=
    Add thermodynamic feasibility constraint (ΔG < 0 for a feasible reaction in flux direction).
    Add objective constraint to solve max min problem.
    =#
    constraints *=
        :reaction_delta_G_margin^C.ConstraintTree(
            Symbol(rxn) => C.Constraint(
                value = constraints.delta_G_reactions[Symbol(rxn)].value *
                        sign(get(reference_flux, rxn, 1.0)),
                bound = C.Between(-Inf, 0.0),
            ) for rxn in rxns
        )

    constraints *=
        :min_driving_force_margin^C.ConstraintTree(
            Symbol(rxn) => C.Constraint(
                value = constraints.max_min_driving_force.value +
                        constraints.delta_G_reactions[Symbol(rxn)].value *
                        sign(get(reference_flux, rxn, 1.0)),
                bound = C.Between(-Inf, 0.0),
            ) for rxn in rxns
        )

    constraints
end

export max_min_driving_force_constraints
