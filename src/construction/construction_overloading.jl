struct MetaboliteWithCoefficient
    coeff :: Float64
    metabolite :: Metabolite
    MetaboliteWithCoefficient(c, m) = new(Float64(c), m)
end

function Base.:*(coeff, met::Metabolite)
    return MetaboliteWithCoefficient(coeff, met)
end

function Base.:+(m1::Union{Metabolite, MetaboliteWithCoefficient}, m2::Union{Metabolite, MetaboliteWithCoefficient})
    if typeof(m1) == Metabolite
        m1 = MetaboliteWithCoefficient(1.0, m1)
    end
    if typeof(m2) == Metabolite
        m2 =  MetaboliteWithCoefficient(1.0, m2)
    end
    return MetaboliteWithCoefficient[m1, m2]
end

function Base.:+(m1::Array{MetaboliteWithCoefficient, 1}, m2::Union{Metabolite, MetaboliteWithCoefficient})
    if typeof(m2) == Metabolite
        m2 = MetaboliteWithCoefficient(1.0, m2)
    end
    return push!(m1, m2)
end

function mkrxn(substrates, products)
    metdict = Dict{Metabolite, Float64}()
    
    if typeof(substrates) == Metabolite
        substrates != ∅ && (metdict[substrates] = get(metdict, substrates, 0.0)  - 1.0) 
    elseif typeof(substrates) == MetaboliteWithCoefficient
        substrates.metabolite != ∅ && (metdict[substrates.metabolite] = get(metdict, substrates.metabolite, 0.0) - 1.0*abs(substrates.coeff))
    else
        for mwc in substrates
            metdict[mwc.metabolite] = get(metdict, mwc.metabolite, 0.0) - 1.0*abs(mwc.coeff)
        end
    end

    if typeof(products) == Metabolite
        products != ∅ && (metdict[products] = get(metdict, products, 0.0) + 1.0) 
    elseif typeof(products) == MetaboliteWithCoefficient
        products.metabolite != ∅ && (metdict[products.metabolite] = get(metdict, products.metabolite, 0.0) + abs(products.coeff))
    else
        for mwc in products
            metdict[mwc.metabolite] = get(metdict, mwc.metabolite, 0.0) + 1.0*abs(mwc.coeff)
        end
    end
    
    return metdict
end

"""
Forward reaction.
"""
function ⟶(substrates::Union{Metabolite, MetaboliteWithCoefficient, Array{MetaboliteWithCoefficient, 1}}, products::Union{Metabolite, MetaboliteWithCoefficient, Array{MetaboliteWithCoefficient, 1}})
    metdict = mkrxn(substrates, products)
    return Reaction("", metdict, "for")
end
const → = ⟶

"""
Reverse only reaction.
"""
function ⟵(substrates::Union{Metabolite, MetaboliteWithCoefficient, Array{MetaboliteWithCoefficient, 1}}, products::Union{Metabolite, MetaboliteWithCoefficient, Array{MetaboliteWithCoefficient, 1}})
    metdict = mkrxn(substrates, products)
    return Reaction("", metdict, "rev")
end
const ← = ⟵

"""
Bidirectional (reversible) reaction.
"""
function ⟷(substrates::Union{Metabolite, MetaboliteWithCoefficient, Array{MetaboliteWithCoefficient, 1}}, products::Union{Metabolite, MetaboliteWithCoefficient, Array{MetaboliteWithCoefficient, 1}})
    metdict = mkrxn(substrates, products)
    return Reaction("", metdict, "bidir")
end
const ↔ = ⟷ 
const ⇌ = ⟷