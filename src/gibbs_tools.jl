"""
Standard Gibbs free energy of reaction
"""
struct ΔᵣG⁰
    ph :: Float64
    rxn :: Union{Reaction, Nothing} # get stoichiometry and concentrations from here
    ΔG :: Union{Float64, Measurement{Float64}} 
end

"""
Adjusted Gibbs free energy of reaction. 

ΔᵣG = ΔᵣG⁰ + RTln(Q)
"""
mutable struct ΔᵣG
    R :: Float64
    T :: Float64
    Q :: Float64 # this gets calculated.
    ΔG⁰ :: Union{Float64, Measurement{Float64}} 
    ΔG :: Union{Float64, Measurement{Float64}} 
end

function update!(ΔG::ΔᵣG)
    if !isnothing(ΔG.ΔG⁰.rxn)
        for (met, stoich) in ΔG.ΔG⁰.rxn
            if met.concentration ≈ 0.0 # no zero concentrations
                ΔG.Q = 1.0            
                break
            else
                ΔG.Q = stoich > 0 ? (ΔG.Q = ΔG.Q*met.concentration^stoich) : (ΔG.Q = ΔG.Q/met.concentration^stoich)
            end
        end
    else
        ΔG.Q = 1.0
    end
    ΔG.ΔG = ΔG.ΔG⁰ + ΔG.R * ΔG.T * log(ΔG.Q)  
end

"""
Display standard ΔG.
"""
function Base.show(io::IO, ΔG::ΔᵣG⁰)
    println(io, "ΔᵣG⁰  (pH = $(ΔG.pH))= ", dg.ΔG ," kJ/mol")
end

"""
Display modified ΔG.
"""
function Base.show(io::IO, ΔG::ΔᵣG)
    println(io, "ΔᵣG⁰  (pH = $(ΔG.pH))= ", dg.ΔG ," kJ/mol")
    println(io, "At T = $(ΔG.T) and Q = $(ΔG.Q)")
    if !isnothing(ΔG.ΔG.rxn)
        for (met, stoich) in ΔG.ΔG.rxn
            println(io, met.name, " = ", met.concentration)
        end
    end
end