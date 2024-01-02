
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

Return true if the reaction denoted by `rxn_id` in `model` is a boundary
reaction, otherwise return false. Checks if on boundary by inspecting the number
of metabolites in the reaction stoichiometry. Boundary reactions have only one
metabolite, e.g. an exchange reaction, or a sink/demand reaction.
"""
is_boundary(model::A.AbstractFBCModel, rxn_id::String) =
    length(keys(A.reaction_stoichiometry(model, rxn_id))) == 1

export is_boundary
