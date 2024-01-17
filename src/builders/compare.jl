
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

TODO
"""
difference_constraint(a, b, distance_bound) =
    C.Constraint(C.value(b) - C.value(a), distance)

"""
$(TYPEDSIGNATURES)

TODO
"""
same_value_constraint(a, b) = C.Constraint(C.value(a) - C.value(b), 0)

"""
$(TYPEDSIGNATURES)

TODO
"""
all_same_constraints(a, bs::C.ConstraintTree) =
    C.map(bs) do b
        same_value_constraint(a, b)
    end

"""
$(TYPEDSIGNATURES)

TODO
"""
greater_or_equal_constraint(a, b) = C.Constraint(C.value(a) - C.value(b), C.Between(0, Inf))

"""
$(TYPEDSIGNATURES)

TODO
"""
less_or_equal_constraint(a, b) = C.Constraint(C.value(b) - C.value(a), C.Between(0, Inf))

# TODO try to use the helper functions everywhere
