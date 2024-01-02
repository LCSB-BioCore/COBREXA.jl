
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
supplied ΔG⁰s in `reaction_standard_gibbs_free_energies`, as defined by Noor, et
al., "Pathway thermodynamics highlights kinetic obstacles in central
metabolism.", PLoS computational biology, 2014.

Optionally, `reference_flux` can be used to set the directions of each reaction
in `model` (all reactions are assumed to proceed forward by default). The
supplied `reference_flux` should be free of internal cycles i.e.
thermodynamically consistent. This optional input is important if a reaction in
`model` normally runs in reverse (negative flux). Note, only the signs are
extracted from this input, so be careful with floating point precision near 0.

The max-min driving force algorithm returns the Gibbs free energy of the
reactions, the concentrations of metabolites and the actual maximum minimum
driving force. The optimization problem solved is:
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

Reactions specified in `ignore_reaction_ids` are internally ignored when
calculating the max-min driving force. Importantly, this should include water
and proton importers.

Since biochemical thermodynamics are assumed, the `proton_ids` and `water_ids`
need to be specified so that they can be ignored in the calculations.
Effectively this assumes an aqueous environment at constant pH is used.

`constant_concentrations` is used to fix the concentrations of certain
metabolites (such as CO₂) by passing a dictionary mapping metabolite id to its
constant value. `concentration_ratios` is used to specify additional constraints
on metabolite pair concentrations (typically, this is done with various
cofactors, such as the ATP/ADP ratio. For example, you can fix the concentration
of ATP to be always 5× higher than of ADP by specifying `Dict("atp_ratio" =>
("atp_c","adp_c", 5.0))` if these metabolites are called `atp_c` and `adp_c` in
the model. `concentration_lb` and `concentration_ub` set the `Cₗ` and `Cᵤ` in
the optimization problems (these are overwritten by `constant_concentrations` if
supplied).

`T` and `R` can be specified in the corresponding units; defaults are K and
kJ/K/mol. The unit of metabolite concentrations is typically molar, and the ΔG⁰s
have units of kJ/mol. Other units can be used, as long as they are consistent.
As usual, optimizer settings can be changed with `modifications`.
"""
function max_min_driving_force_analysis(
    model::A.AbstractFBCModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64};
    reference_flux = Dict{String,Float64}(),
    concentration_ratios = Dict{String,Tuple{String,String,Float64}}(),
    constant_concentrations = Dict{String,Float64}(),
    proton_ids,
    water_ids,
    concentration_lb = 1e-9, # M
    concentration_ub = 1e-1, # M
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = String[],
    modifications = [],
    optimizer,
)
    m = build_max_min_driving_force_model(
        model,
        reaction_standard_gibbs_free_energies;
        reference_flux,
        concentration_lb,
        concentration_ub,
        R,
        T,
        ignore_reaction_ids,
        water_ids,
        proton_ids,
    )

    for (mid, val) in constant_concentrations
        m.log_metabolite_concentrations[Symbol(mid)].bound = C.EqualTo(log(val))
    end

    m *=
        :metabolite_ratio_constraints^log_ratio_constraints(
            concentration_ratios,
            m.log_metabolite_concentrations,
        )


    optimized_constraints(
        m;
        objective = m.max_min_driving_force.value,
        optimizer,
        modifications,
    )
end

export max_min_driving_force_analysis
