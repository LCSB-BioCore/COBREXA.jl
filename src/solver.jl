
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
$(TYPEDEF)

Representation of a "binary switch" bound for `ConstraintTree`s. The value is
constrained to be either the value of field `a` or of field `b`; both fields
are `Float64`s. Upon translation to JuMP, the switches create an extra boolean
variable, and the value is constrained to equal `a + boolean_var * (b-a)`.

Switches can be offset by adding real numbers, negated, and multiplied and
divided by scalar constraints. For optimizing some special cases, multiplying
by exact zero returns an equality bound to zero.
"""
struct Switch <: C.Bound
    a::Float64
    b::Float64
end

export Switch

Base.:-(x::Switch) = Switch(-s.a, -s.b)
Base.:+(x::Real, s::Switch) = b + a
Base.:+(s::Switch, x::Real) = Switch(s.a + x, s.b + x)
Base.:*(x::Real, s::Switch) = b * a
Base.:*(s::Switch, x::Real) = x == 0 ? C.EqualTo(0) : Switch(s.a * x, s.b * x)
Base.:/(s::Switch, x::Real) = Switch(s.a / x, s.b / x)

"""
$(TYPEDSIGNATURES)

Construct a JuMP `Model` that describes the precise constraint system into the
JuMP `Model` created for solving in `optimizer`, with a given optional
`objective` and optimization `sense` chosen from [`Maximal`](@ref),
[`Minimal`](@ref) and [`Feasible`](@ref).
"""
function optimization_model(
    cs::C.ConstraintTreeElem;
    objective::Union{Nothing,C.Value} = nothing,
    optimizer,
    sense = Maximal,
)
    model = J.Model(optimizer)

    J.@variable(model, x[1:C.var_count(cs)])
    isnothing(objective) || J.@objective(model, sense, C.substitute(objective, x))

    # constraints
    function add_constraint(v::C.Value, b::C.EqualTo)
        J.@constraint(model, C.substitute(v, x) == b.equal_to)
    end
    function add_constraint(v::C.Value, b::C.Between)
        vx = C.substitute(v, x)
        isinf(b.lower) || J.@constraint(model, vx >= b.lower)
        isinf(b.upper) || J.@constraint(model, vx <= b.upper)
    end
    function add_constraint(v::C.Value, b::Switch)
        boolean = J.@variable(model, binary = true)
        J.@constraint(model, C.substitute(v, x) == b.a + boolean * (b.b - b.a))
    end
    add_constraint(::C.Value, _::Nothing) = nothing
    function add_constraint(c::C.Constraint)
        add_constraint(c.value, c.bound)
    end
    function add_constraint(c::C.ConstraintTree)
        add_constraint.(values(c))
    end

    add_constraint(cs)

    return model
end

export optimization_model

"""
$(TYPEDSIGNATURES)

`true` if `opt_model` solved successfully (solution is optimal or
locally optimal). `false` if any other termination status is reached.
"""
is_solved(opt_model::J.Model) =
    J.termination_status(opt_model) in [J.MOI.OPTIMAL, J.MOI.LOCALLY_SOLVED]

export is_solved

"""
    Minimal

Objective sense for finding the minimal value of the objective.

Same as `JuMP.MIN_SENSE`.
"""
const Minimal = J.MIN_SENSE
export Minimal

"""
    Maximal

Objective sense for finding the maximal value of the objective.

Same as `JuMP.MAX_SENSE`.
"""
const Maximal = J.MAX_SENSE
export Maximal

"""
    Maximal

Objective sense for finding the any feasible value of the objective.

Same as `JuMP.FEASIBILITY_SENSE`.
"""
const Feasible = J.FEASIBILITY_SENSE
export Feasible

"""
$(TYPEDSIGNATURES)

Make an JuMP model out of `constraints` using [`optimization_model`](@ref)
(most arguments are forwarded there), then apply the `settings`, optimize
the model, and return either `nothing` if the optimization failed, or `output`
substituted with the solved values (`output` defaults to `constraints`.

For a "nice" version for simpler finding of metabolic model optima, use
[`flux_balance`](@ref).
"""
function optimized_constraints(
    constraints::C.ConstraintTreeElem;
    settings = [],
    output = constraints,
    kwargs...,
)
    om = optimization_model(constraints; kwargs...)
    for m in settings
        m(om)
    end
    J.optimize!(om)
    is_solved(om) ? C.substitute_values(output, J.value.(om[:x])) : nothing
end

export optimized_constraints
