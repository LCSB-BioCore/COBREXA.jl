struct Gibbs 
    ΔG :: Measurement{Float64}
end

mutable struct ReactionGibbs
    rxns :: Array{Reaction, 1}    
    ΔG⁰s :: Array{Gibbs, 1} # Gibbs data - pH dependent
    ΔGs :: Array{Gibbs, 1} # Adjusted by concentration, temperature
end

function Base.show(io::IO, rgs::ReactionGibbs)
    println(io, "Gibbs free energy of reaction set of length: $(length(rgs.rxns))")
    println(io, "Fraction non-zero: $(count(x->x != 0.0, [x.ΔG for x in rgs.ΔG⁰s])/length(rgs.rxns))")
     
end

function Base.show(io::IO, gibbs::Gibbs)
    println(io, "ΔᵣG⁰ = ", gibbs.ΔG, " kJ/mol")
end

function mapGibbs(rxns::Array{Reaction, 1}; ph="7.0", T=298.15, R=8314.46261815324)
    CC = pyimport("from equilibrator_api import ComponentContribution")
    Q = pyimport("from equilibrator_api import Q_")
    # cc = ComponentContribution()
end

function getdG0(formula::String, ph=7.0, ionic_strength="100 mM")
    isbal, v, err = py"pygetdg0"(formula, ph, ionic_strength)
    !isbal && @warn "Reaction not balanced."
    return v ± err
end