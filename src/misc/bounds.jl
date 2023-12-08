
"""
$(TYPEDSIGNATURES)

TODO
"""
absolute_tolerance_bound(tolerance) = x -> begin
    bound = (x - tolerance, x + tolerance)
    (minimum(bound), maximum(bound))
end

"""
$(TYPEDSIGNATURES)

TODO
"""
relative_tolerance_bound(tolerance) = x -> begin
    bound = (x * tolerance, x / tolerance)
    (minimum(bound), maximum(bound))
end
