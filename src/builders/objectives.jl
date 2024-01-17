
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

Construct a `ConstraintTrees.Value` out of squared sum of all values directly
present in a given constraint tree.
"""
squared_sum_value(x::C.ConstraintTree) = squared_sum_error_value(x, Dict(keys(x) .=> 0.0))

"""
$(TYPEDSIGNATURES)

Construct a `ConstraintTrees.Value` out of a sum of all values directly present
in a given constraint tree.
"""
function sum_value(x...)
    res = zero(C.LinearValue)
    for ct in x
        C.map(ct) do c
            res += c.value
        end
    end
    res
end

"""
$(TYPEDSIGNATURES)

Construct a `ConstraintTrees.Value` out of squared error (in the RMSE-like
squared-error sense) between the values in the constraint tree and the
reference `target`.

`target` is a function that takes a symbol (key) and returns either a `Float64`
reference value, or `nothing` if the error of given key should not be
considered.
"""
squared_sum_error_value(constraints::C.ConstraintTree, target) = sum(
    (
        C.squared(C.value(c) - t) for
        (t, c) in ((target(k), c) for (k, c) in constraints) if !isnothing(t)
    ),
    init = zero(C.LinearValue),
)
