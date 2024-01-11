
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

Add loopless constraints to the model, `m`. Specify the internal reactions with
`internal_reaction_ids`, as well as the columns of the stoichiometric nullspace
of these reactions in `internal_reaction_stoichiometry_nullspace_columns`. See
the example below for more information about the nullspace. By default, the
`fluxes` of the model `m` are made loopless.

The Big-M method is used to ensure the sign of fluxes and pseudo Gibbs free
energy of reactions match. For this, ensure that `max_flux_bound` is at least
one order of magnitude bigger than the largest expected absolute flux value.
Additionally, ensure that `strict_inequality_tolerance` is smaller than any
expected pseudo Gibbs free energy of reaction value. The defaults work well for
most problems.

# Example
```
internal_reaction_stoichiometry_nullspace_columns =
    eachcol(nullspace(Array(A.stoichiometry(model)[:, internal_rxn_idxs_in_order_of_internal_rxn_ids])))
```
"""
function loopless_constraints(;
    fluxes::C.ConstraintTree,
    loopless_direction_indicators::C.ConstraintTree,
    loopless_driving_forces::C.ConstraintTree,
    internal_reactions::Vector{Symbol},
    internal_nullspace::Matrix,
    flux_infinity_bound,
    driving_force_nonzero_bound,
    driving_force_infinity_bound,
)

    C.ConstraintTree(
        :flux_direction_lower_bounds => C.ConstraintTree(
            r => C.Constraint(
                value = fluxes[r].value +
                        flux_infinity_bound * (1 - loopless_direction_indicators[r].value),
                bound = C.Between(0, Inf),
            ) for r in internal_reactions
        ),
        :flux_direction_upper_bounds => C.ConstraintTree(
            r => C.Constraint(
                value = fluxes[r].value +
                        flux_infinity_bound * loopless_direction_indicators[r].value,
                bound = C.Between(-Inf, 0),
            ) for r in internal_reactions
        ),
        :driving_force_lower_bounds => C.ConstraintTree(
            r => C.Constraint(
                value = loopless_driving_forces[r].value -
                        strict_inequality_tolerance *
                        loopless_direction_indicators[r].value +
                        flux_infinity_bound * (1 - loopless_direction_indicators[r].value),
                bound = C.Between(0, Inf),
            ) for r in internal_reaction_ids
        ),
        :driving_force_upper_bounds => C.ConstraintTree(
            r => C.Constraint(
                value = loopless_driving_forces[r].value +
                        strict_inequality_tolerance *
                        (1 - loopless_direction_indicators[r].value) -
                        flux_infinity_bound * loopless_direction_indicators[r].value,
                bound = C.Between(-Inf, 0),
            ) for r in internal_reaction_ids
        ),
        :loopless_nullspace => C.ConstraintTree(
            Symbol(:nullspace_base_, i) => C.Constraint(
                value = sum(
                    coeff * loopless_driving_forces[r].value for
                    (coeff, r) in zip(vec, internal_reactions)
                ),
                bound = C.EqualTo(0),
            ) for (i, col) in enumerate(eachcol(internal_nullspace))
        ),
    )
end

export loopless_constraints
