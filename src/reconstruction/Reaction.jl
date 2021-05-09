struct MetaboliteWithCoefficient
    coeff::Float64
    metabolite::Metabolite
    MetaboliteWithCoefficient(c, m) = new(Float64(c), m)
end

function Base.:*(coeff, met::Metabolite)
    return MetaboliteWithCoefficient(coeff, met)
end

function Base.:+(
    m1::Union{Metabolite,MetaboliteWithCoefficient},
    m2::Union{Metabolite,MetaboliteWithCoefficient},
)
    if typeof(m1) == Metabolite
        m1 = MetaboliteWithCoefficient(1.0, m1)
    end
    if typeof(m2) == Metabolite
        m2 = MetaboliteWithCoefficient(1.0, m2)
    end
    return MetaboliteWithCoefficient[m1, m2]
end

function Base.:+(
    m1::Vector{MetaboliteWithCoefficient},
    m2::Union{Metabolite,MetaboliteWithCoefficient},
)
    if typeof(m2) == Metabolite
        m2 = MetaboliteWithCoefficient(1.0, m2)
    end
    return push!(m1, m2)
end

function _mkrxn(substrates, products)
    metdict = Dict{String,Float64}()

    if typeof(substrates) == Metabolite
        metdict[substrates.id] = get(metdict, substrates.id, 0.0) - 1.0
    elseif typeof(substrates) == MetaboliteWithCoefficient
        metdict[substrates.metabolite.id] =
            get(metdict, substrates.metabolite.id, 0.0) - 1.0 * abs(substrates.coeff)
    elseif typeof(products) == Vector{MetaboliteWithCoefficient}
        for mwc in substrates
            metdict[mwc.metabolite.id] =
                get(metdict, mwc.metabolite.id, 0.0) - 1.0 * abs(mwc.coeff)
        end
    end

    if typeof(products) == Metabolite
        metdict[products.id] = get(metdict, products.id, 0.0) + 1.0
    elseif typeof(products) == MetaboliteWithCoefficient
        metdict[products.metabolite.id] =
            get(metdict, products.metabolite.id, 0.0) + abs(products.coeff)
    elseif typeof(products) == Vector{MetaboliteWithCoefficient}
        for mwc in products
            metdict[mwc.metabolite.id] =
                get(metdict, mwc.metabolite.id, 0.0) + 1.0 * abs(mwc.coeff)
        end
    end

    return metdict
end

"""
    ⟶(
        substrates::Union{
            Nothing,
            Metabolite,
            MetaboliteWithCoefficient,
            Vector{MetaboliteWithCoefficient},
        },
        products::Union{
            Nothing, 
            Metabolite,
            MetaboliteWithCoefficient,
            Vector{MetaboliteWithCoefficient}
        },
    )

Make a forward-only [`Reaction`](@ref) from `substrates` and `products`.
An equivalent alternative is `→`.
"""
function ⟶(
    substrates::Union{
        Nothing,
        Metabolite,
        MetaboliteWithCoefficient,
        Vector{MetaboliteWithCoefficient},
    },
    products::Union{
        Nothing, 
        Metabolite,
        MetaboliteWithCoefficient,
        Vector{MetaboliteWithCoefficient}
    },
)
    metdict = _mkrxn(substrates, products)
    return Reaction("", metdict, :forward)
end
const → = ⟶

"""
    ⟵(
        substrates::Union{
            Nothing,
            Metabolite,
            MetaboliteWithCoefficient,
            Vector{MetaboliteWithCoefficient},
        },
        products::Union{
            Nothing,
            Metabolite,
            MetaboliteWithCoefficient,
            Vector{MetaboliteWithCoefficient}
        },
    )

Make a reverse-only [`Reaction`](@ref) from `substrates` and `products`.
An equivalent alternative is `←`.
"""
function ⟵(
    substrates::Union{
        Nothing,
        Metabolite,
        MetaboliteWithCoefficient,
        Vector{MetaboliteWithCoefficient},
    },
    products::Union{
        Nothing,
        Metabolite,
        MetaboliteWithCoefficient,
        Vector{MetaboliteWithCoefficient}
    },
)
    metdict = _mkrxn(substrates, products)
    return Reaction("", metdict, :reverse)
end
const ← = ⟵

"""
    ⟷(
        substrates::Union{
            Nothing,
            Metabolite,
            MetaboliteWithCoefficient,
            Vector{MetaboliteWithCoefficient},
        },
        products::Union{
            Nothing,
            Metabolite,
            MetaboliteWithCoefficient,
            Vector{MetaboliteWithCoefficient}
        },
    )

Make a bidirectional (reversible) [`Reaction`](@ref) from `substrates` and `products`.
An equivalent alternative is `↔`.
"""
function ⟷(
    substrates::Union{
        Nothing,
        Metabolite,
        MetaboliteWithCoefficient,
        Vector{MetaboliteWithCoefficient},
    },
    products::Union{
        Nothing,
        Metabolite,
        MetaboliteWithCoefficient,
        Vector{MetaboliteWithCoefficient}
    },
)
    metdict = _mkrxn(substrates, products)
    return Reaction("", metdict, :bidirectional)
end
const ↔ = ⟷
