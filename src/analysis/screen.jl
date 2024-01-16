
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

function screen(
    f,
    model::A.AbstractFBCModel;
    args::Array{Tuple} = [()],
    workers = D.workers(),
    kwargs...,
)
    # TODO might belong to the frontend
end

function screen(
    f,
    constraints::C.ConstraintTree;
    args::Array{Tuple} = [()],
    workers = D.workers(),
)
    # TODO
end

function screen_optimization_model(
    f,
    constraints::C.ConstraintTree,
    args::Maybe{Array},
    workers = D.workers(),
)
    # TODO
end
