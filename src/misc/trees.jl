
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

#TODO these are likely hot candidates to be moved to CTs

"""
$(TYPEDSIGNATURES)

Extract all elements of a `ConstraintTrees.Tree` in order and return them in a
`Vector` transformed by `f`. If the order is not modified, one can re-insert a
vector of modified elements into the same-shaped tree using
[`tree_reinflate`](@ref).
"""
function tree_deflate(f, x::C.Tree{T})::Vector{T} where {T}
    count = 0
    C.traverse(x) do _
        count += 1
    end
    res = Vector{T}(undef, count)
    i = 1
    C.traverse(x) do c
        res[i] = f(c)
        i += 1
    end
    res
end

"""
$(TYPEDSIGNATURES)

Insert a `Vector` of elements into the "values" of a `ConstraintTrees.Tree`.
The order of elements is given by [`tree_deflate`](@ref).
"""
function tree_reinflate(x::C.Tree, elems::Vector{T})::C.Tree{T} where {T}
    i = 0
    C.map(x) do _
        i += 1
        elems[i]
    end
end
