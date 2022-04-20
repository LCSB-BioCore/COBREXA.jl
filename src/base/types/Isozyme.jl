"""
    mutable struct Isozyme

Information about isozyme composition and activity.

# Fields
- `gene_product_count :: Dict{String, Int}` assigns each gene product ID its
  count in the isozyme complex (which is used to determine the total mass of
  the isozyme)
- `kcat_forward`, `kcat_reverse` -- forward and reverse turnover numbers of the
  isozyme
````
"""
mutable struct Isozyme
    gene_product_count::Dict{String,Int}
    kcat_forward::Float64
    kcat_reverse::Float64
end
