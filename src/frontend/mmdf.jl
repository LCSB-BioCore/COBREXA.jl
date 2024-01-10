
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

Perform a max-min driving force analysis using `optimizer` on the `model` with
supplied reaction standard Gibbs energies in `reaction_standard_gibbs_free_energies`.

The method was described by by: Noor, et al., "Pathway thermodynamics highlights
kinetic obstacles in central metabolism.", PLoS computational biology, 2014.

`reference_flux` sets the directions of each reaction in `model`. The scale of
the values is not important, only the direction is examined (w.r.t.
`reference_flux_atol` tolerance). Ideally, the supplied `reference_flux` should
be completely free of internal cycles, which enables the thermodynamic
consistency. To get the cycle-free flux, you can use
[`loopless_flux_balance_analysis`](@ref) (computationally complex but gives
very good fluxes), [`parsimonious_flux_balance_analysis`](@ref) or
[`linear_parsimonious_flux_balance_analysis`](@ref) (computationally simplest
but the consistency is not guaranteed).

Internally, [`fbc_log_concentration_constraints`](@ref) is used to lay out the
base structure of the problem.

Following arguments are set optionally:
- `water_metabolites`, `proton_metabolites` and `ignored_metabolites` allow to
  completely ignore constraints on a part of the metabolite set, which is
  explicitly recommended especially for water and protons (since the analyses
  generally assume aqueous environment of constant pH)
- `constant_concentrations` can be used to fix the concentrations of the
  metabolites
- `concentration_lower_bound` and `concentration_upper_bound` set the default
  concentration bounds for all other metabolites
- `concentration ratios` is a dictionary that assigns a tuple of
  metabolite-metabolite-concentration ratio constraint names; e.g. ATP/ADP
  ratio can be fixed to five-times-more-ATP by setting `concentration_ratios =
  Dict("adenosin_ratio" => ("atp", "adp", 5.0))`
- `T` and `R` default to the usual required thermodynamic constraints in the
  usual units (K and kJ/K/mol, respectively). These multiply the
  log-concentrations to obtain the actual Gibbs energies and thus driving
  forces.
"""
function max_min_driving_force_analysis(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    reference_flux = Dict{String,Float64}(),
    concentration_ratios = Dict{String,Tuple{String,String,Float64}}(),
    constant_concentrations = Dict{String,Float64}(),
    ignored_metabolites = Set{String}(),
    proton_metabolites = Set{String}(),
    water_metabolites = Set{String}(),
    concentration_lower_bound = 1e-9, # M
    concentration_upper_bound = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    reference_flux_atol = 1e-6,
    check_ignored_reactions = missing,
    settings = [],
    optimizer,
)

    # First let's just check if all the identifiers are okay because we use
    # quite a lot of these; the rest of the function may be a bit cleaner with
    # this checked properly.

    model_reactions = Set(A.reactions(model))
    model_metabolites = Set(A.metabolites(model))

    all(in(model_reactions), keys(reaction_standard_gibbs_free_energies)) || throw(
        DomainError(
            reaction_standard_gibbs_free_energies,
            "unknown reactions referenced by reaction_standard_gibbs_free_energies",
        ),
    )
    all(x -> haskey(reaction_standard_gibbs_free_energies, x), keys(reference_flux)) ||
        throw(DomainError(reference_flux, "some reactions have no reference flux"))
    all(in(model_reactions), keys(reference_flux)) || throw(
        DomainError(
            reaction_standard_gibbs_free_energies,
            "unknown reactions referenced by reference_flux",
        ),
    )
    all(in(model_metabolites), keys(constant_concentrations)) || throw(
        DomainError(
            constant_concentrations,
            "unknown metabolites referenced by constant_concentrations",
        ),
    )
    all(
        in(model_metabolites),
        (m for (_, (x, y, _)) in concentration_ratios for m in (x, y)),
    ) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by concentration_ratios",
        ),
    )
    all(in(model_metabolites), proton_metabolites) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by proton_metabolites",
        ),
    )
    all(in(model_metabolites), water_metabolites) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by water_metabolites",
        ),
    )
    all(in(model_metabolites), ignored_metabolites) || throw(
        DomainError(
            concentration_ratios,
            "unknown metabolites referenced by ignored_metabolites",
        ),
    )
    if !ismissing(check_ignored_reactions) && (
        all(
            x -> !haskey(reaction_standard_gibbs_free_energies, x),
            check_ignored_reactions,
        ) || (
            union(
                Set(check_ignored_reactions),
                Set(keys(reaction_standard_gibbs_free_energies)),
            ) != model_reactions
        )
    )
        throw(AssertionError("check_ignored_reactions validation failed"))
    end

    # that was a lot of checking.

    default_concentration_bound =
        C.Between(log(concentration_lower_bound), log(concentration_upper_bound))

    no_concentration_metabolites = union(
        Set(Symbol.(water_metabolites)),
        Set(Symbol.(proton_metabolites)),
        Set(Symbol.(ignored_metabolites)),
    )

    constraints =
        fbc_log_concentration_constraints(
            model,
            concentration_bound = met -> if met in no_concentration_metabolites
                C.EqualTo(0)
            else
                mid = String(met)
                if haskey(constant_concentrations, mid)
                    C.EqualTo(log(constant_concentrations[mid]))
                else
                    default_concentration_bound
                end
            end,
        ) + :max_min_driving_force^C.variable()

    min_driving_forces = C.ConstraintTree(
        let r = Symbol(rid)
            r => C.Constraint(
                value = dGr0 + (R * T) * constraints.reactant_log_concentrations[r],
                bound = let rf = reference_flux[rid]
                    if isapprox(rf, 0.0, atol = reference_flux_atol)
                        C.EqualTo(0)
                    elseif rf > 0
                        C.Between(-Inf, 0)
                    else
                        C.Between(0, Inf)
                    end
                end,
            )
        end for (rid, dGr0) in reaction_standard_gibbs_free_energies
    )

    constraints =
        constraints *
        :min_driving_forces^min_driving_forces *
        :min_driving_force_thresholds^C.map(min_driving_forces) do c
            C.Constraint(
                value = constraints.max_min_driving_force.value - c.value;
                bound = C.Between(0, Inf),
            )
        end *
        :concentration_ratio_constraints^C.ConstraintTree(
            Symbol(cid) => C.Constraint(
                m.log_metabolite_concentrations[Symbol(m2)].value -
                m.log_metabolite_concentrations[Symbol(m1)].value,
                log(ratio),
            ) for (cid, (m1, m2, ratio)) in concentration_ratios
        )

    optimized_constraints(m; objective = m.max_min_driving_force.value, optimizer, settings)
end

export max_min_driving_force_analysis
