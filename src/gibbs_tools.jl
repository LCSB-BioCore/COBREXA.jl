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
    gibbs_data = JLD.load(joinpath("data", "gibbs_$ph.jld"), "gibbs")
    name_mapping = JLD.load(joinpath("data", "bigg_to_kegg.jld"), "bigg_to_kegg") # need to fix for other name spaces
    ΔG⁰s = Gibbs[]
    ΔGs = Gibbs[]
    for i in eachindex(rxns)
        id = rxns[i].id
        if haskey(name_mapping, id)
            gs = Measurement{Float64}[]
            for kegg in name_mapping[id]
                if haskey(gibbs_data, kegg)
                    push!(gs, gibbs_data[kegg][1] ± gibbs_data[kegg][2])
               end
            end
            if isempty(gs)
                println(id)
                push!(ΔG⁰s, Gibbs(0.0 ± 0.0)) # missing 
                push!(ΔGs, Gibbs(0.0 ± 0.0)) # missing    
            else
                mg = mean(gs)
                Q = 1.0
                for (met, stoich) in rxns[1].metabolites
                    if stoich > 0
                        Q = Q * (met.concentration^stoich)
                    else
                        Q = Q / (met.concentration^stoich)
                    end
                end
                gadj = mg + R * T * log(Q)
                push!(ΔG⁰s, Gibbs(mg)) # missing 
                push!(ΔGs, Gibbs(gadj)) # missing    
            end
        else
            println(id)
            push!(ΔG⁰s, Gibbs(0.0 ± 0.0)) # missing 
            push!(ΔGs, Gibbs(0.0 ± 0.0)) # missing
        end
    end
    return ReactionGibbs(rxns, ΔG⁰s, ΔGs)
end