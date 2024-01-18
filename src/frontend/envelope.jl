
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

Find the objective production envelope of the `model` in the dimensions given
by `reactions`.

This runs a variability analysis of constraints to determine an applicable
range for the dimensions, then splits the dimensions to equal-sized breaks (of
count `breaks` for each dimension, i.e. total `breaks ^ length(reactions)`
individual "multidimensional breaks") thus forming a grid, and returns an array
of fluxes through the model objective with the individual reactions fixed to
flux as given by the grid.

`optimizer` and `settings` are used to construct the optimization models.

The computation is distributed to the specified `workers`, defaulting to all
available workers.
"""
function objective_production_envelope(
    model::A.AbstractFBCModel,
    reactions::Vector{String};
    breaks = 10,
    optimizer,
    settings = [],
    workers = D.workers(),
)
    constraints = flux_balance_constraints(model)
    rs = Symbol.(reactions)

    envelope_bounds = constraints_variability(
        constraints,
        (r => constraints.fluxes[r] for r in rs);
        optimizer,
        settings,
        workers,
    )

    #TODO check for nothings in the bounds

    bss = [break_interval(envelope_bounds[r]..., breaks) for r in rs]

    return (
        breaks = reactions .=> bss,
        objective_values = constraints_objective_envelope(
            constraints,
            (constraints.fluxes[r] => bs for (r, bs) in zip(rs, bss))...;
            objective = model.objective.value,
            sense = Maximal,
            optimizer,
            settings,
            workers,
        ),
    )

    # this converts nicely to a dataframe, but I'm not a total fan.
    #=
    xss = Iterators.product(bss)
    @assert length(result) == length(xss)
    xss = reshape(xss, tuple(length(xss)))

    return (;
        (r => [xs[i] for xs in xss] for (i, r) in enumerate(rs))...,
        objective_value_name => reshape(result, tuple(length(result))),
    )
    # TODO think about it
    =#
end
