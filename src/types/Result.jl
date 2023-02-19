"""
$(TYPEDEF)

A simple storage type for analysis result that is accompanied with the original
model. This vastly simplifies piping; e.g., instead of having to specify the
model in each part of the pipeline as such:
```
flux_balance_analysis(m, optimizer) |> values_vec(m)
```
...we can do:
```
flux_balance_analysis(m, optimizer) |> values_vec
```

Optionally you may take out "just" the result value by piping through
[`result`](@ref), in this case obtaining a vector:
```
... |> values_vec |> result
```

This additionally enables some situations where carrying the model around
manually would be hard or require multiple piping steps, such as:
```
model |> with_some_reactions(...) |>
         with_enzyme_constraints(...) |>
         flux_balance_analysis(optimizer) |>
         values_dict |>
         result
"""
struct ModelWithResult{T}
    model::AbstractMetabolicModel
    result::T
end

"""
$(TYPEDSIGNATURES)

Pipeable shortcut for extracting the result value out of
[`ModelWithResult`](@ref).
"""
function result(x::ModelWithResult{T})::T where T
    x.result
end
