using CobraTools
using Measurements
using Statistics
using PyCall
using StringDistances
using Plots
gr()

modelpath = joinpath("models", "yeastGEM.mat")
modelpath = joinpath("models", "iJO1366.mat")
modelpath = joinpath("models", "iJO1366.json")

model = CobraTools.read_model(modelpath)
gibbs, balances = CobraTools.map_gibbs_rxns(model.rxns) 
brenda_data = CobraTools.parse_brenda(joinpath("data", "brenda_download.brenda.txt"))

kcats = Dict{String, Union{Float64, Measurement{Float64}}}()
for rxnid in keys(gibbs)
    # rxnid = "GHBDHx"

    metnames = String[]
    rxn = findfirst(model.rxns, rxnid)
    ec = get(rxn.annotation, "ec-code", [""])[1]
    metnames = [m.name for m in keys(rxn.metabolites)]
    br_ind = findfirst(x -> x.ID == ec, brenda_data)
    isnothing(br_ind) && continue

    ks = Float64[]
    for br in brenda_data[br_ind].TN
        isnothing(br.substrate) && continue
        scomps = [compare(lowercase(metname), lowercase(br.substrate), Levenshtein()) for metname in metnames]
        if maximum(scomps) > 0.65 # could probably make this stricter
            push!(ks, br.val)
            # ind = argmax(scomps)
            # println(br[2], " <-> ", metnames[ind], " = ", maximum(scomps))
        end
    end
    if !isempty(ks)
        if length(ks) == 1
            kcats[rxnid] = ks[1]
        else
            kcats[rxnid] = mean(ks) ± std(ks)
        end
    end
end

dkeys = intersect(keys(gibbs), keys(kcats))
plot(xlabel="|ΔG'⁰|", ylabel="Turnover [s⁻¹]")
pvs = []
for dk in dkeys
    if gibbs[dk] == 0.0 # skip default values
        continue 
    end
    ΔG = log10(abs(gibbs[dk]))
    kcat = log10(kcats[dk])
    
    push!(pvs, [ΔG, kcat])
    if typeof(kcat) == Float64
        scatter!([ΔG], [kcat], color="blue")
    else
        if abs(ΔG.err) < 5.0 && abs(kcat.err) < 5.0 
            scatter!([ΔG], [kcat], color="blue")
        end
    end
end
plot!(legend=false) 

plot!(xlim=(-2, 3))
savefig("brenda-dg-kcat.png")