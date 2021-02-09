struct Gibbs 
    ΔG :: Measurement
end

mutable struct ReactionGibbs
    rxns :: Array{Reaction, 1}    
    ΔG⁰s :: Array{Gibbs, 1} # Gibbs data - pH dependent
    ΔGs :: Array{Gibbs, 1} # Adjusted by concentration, temperature
end

function Base.show(io::IO, rgs::ReactionGibbs)
    println(io, "Gibbs free energy of reaction set of length: $(length(rgs))") 
end

function Base.show(io::IO, gibbs::Gibbs)
    println(io, "ΔᵣG⁰ = ", gibbs.ΔG, " kJ/mol")
end

function mapGibbs(rxns::Array{Reaction, 1}; ph="7.0")
    
end