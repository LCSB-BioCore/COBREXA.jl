
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
squared_sum_value(x::C.ConstraintTree) = squared_sum_error_value(x, Dict(keys(x) .=> 0.0))

"""
$(TYPEDSIGNATURES)

TODO useful for L1 parsimonious stuff
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

TODO
"""
squared_sum_error_value(constraints::C.ConstraintTree, target::Dict{Symbol,Float64}) = sum(
    (C.squared(C.value(c) - target[k]) for (k, c) in constraints if haskey(target, k)),
    init = zero(C.LinearValue),
)
