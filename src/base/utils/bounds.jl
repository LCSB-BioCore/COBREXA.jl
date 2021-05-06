"""
    gamma_bounds(gamma)

A bounds-generating function for [`flux_variability_analysis`](@ref) that
limits the objective value to be at least `gamma*Z₀`, as usual in COBRA
packages. Use as the `bounds` argument:
```
flux_variability_analysis(model, some_optimizer; bounds = gamma_bounds(0.9))
```
"""
gamma_bounds(gamma) = z -> (gamma * z, Inf)

"""
    (tolerance) = z -> begin

A bounds-generating function for [`flux_variability_analysis`](@ref) that
limits the objective value to a small multiple of Z₀. Use as `bounds` argument,
similarly to [`gamma_bounds`](@ref).
"""
objective_bounds(tolerance) = z -> begin
    vs = (z * tolerance, z / tolerance)
    (minimum(vs), maximum(vs))
end
