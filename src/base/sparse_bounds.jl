"""
    struct SparseBoundVector
        
A struct used to store bound vectors.

This saves space and is useful in some analysis functions, e.g. sampling.
"""
struct SparseBoundVector{Tv, Ti<:Integer} 
    n::Int # total size of the array
    boundval::Tv # the bound value that occurs the most except for 0
    nonbound_inds::Vector{Ti}
    nonbound_vals::Vector{Tv}
    zero_inds::BitVector
    bound_inds::BitVector # could possibly get rid if this one but it makes the code more complex
end

"""
    sparse_bounds(bound_arr, bound_val)

Return a sparse representation of `bound_arr` where all instances of
`bound_val` and 0 have been removed from `bound_arr`.
"""
function sparse_bounds(arr, bound)
    n = length(arr)
    z_inds = BitVector(undef, n)
    b_inds = BitVector(undef, n)
    v_inds = Int64[]
    vs = Float64[]
    for (i, a) in enumerate(arr)
        if a == 0
            z_inds[i] = 1
            b_inds[i] = 0
        elseif a == bound
            z_inds[i] = 0
            b_inds[i] = 1
        else
            z_inds[i] = 0
            b_inds[i] = 0
            push!(v_inds, i)
            push!(vs, a)
        end
    end
    return SparseBoundVector{Float64, Int64}(n, bound, v_inds, vs, z_inds, b_inds)
end

# Interface is very restricted since the data type is very unique
Base.length(s::SparseBoundVector) = s.n
Base.eltype(s::SparseBoundVector) = typeof(s.nonbound_vals[1])
Base.size(s::SparseBoundVector) = (length(s),)
Base.size(s::SparseBoundVector, d::Int) = length(s)
