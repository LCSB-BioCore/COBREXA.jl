
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

A constraint that makes sure that the difference from `a` to `b` is within the
`difference_bound`. For example, `difference_constraint(-1, 1, 2)` will always
be valid. Any type of `ConstraintTree.Bound` can be supplied.
"""
difference_constraint(a, b, difference_bound) =
    C.Constraint(C.value(b) - C.value(a), difference_bound)

export difference_constraint

"""
$(TYPEDSIGNATURES)

A constraint that makes sure that the values of `a` and `b` are the same.
"""
equal_value_constraint(a, b) = difference_constraint(a, b, 0)

export equal_value_constraint

"""
$(TYPEDSIGNATURES)

A constriant tree that makes sure that all values in `tree` are the same as the
value of `a`.

Names in the output `ConstraintTree` match the names in the `tree`.
"""
all_equal_constraints(a, tree::C.ConstraintTree) =
    C.map(tree) do b
        equal_value_constraint(a, b)
    end

export all_equal_constraints

"""
$(TYPEDSIGNATURES)

A constraint that makes sure that the value of `a` is greater than or equal to
the the value of `b`.
"""
greater_or_equal_constraint(a, b) = difference_constraint(a, b, C.Between(-Inf, 0))

export greater_or_equal_constraint

"""
$(TYPEDSIGNATURES)

A constraint that makes sure that the value of `a` is less than or equal to the
the value of `b`.
"""
less_or_equal_constraint(a, b) = difference_constraint(a, b, C.Between(0, Inf))

export less_or_equal_constraint
