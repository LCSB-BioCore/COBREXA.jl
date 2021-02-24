using CobraTools
using Plots
gr()

# iML1515.json = E. coli K12
# iJN678.json = Synechocystis sp. PCC 6803
# iJN746.json = Pseudomonas putida KT2440
# iRC1080.json = Chlamydomonas reinhardtii
# iYO844.json = Bacillus subtilis subsp. subtilis str. 168
# iMM904.json = S. cerevisiae S288C
# iLJ478.json = Thermotoga maritima MSB8


models = ["iML1515.json", "iJN678.json", "iJN746.json", "iRC1080.json", "iYO844.json", "iMM904.json", "iLJ478.json"]
mnames = Dict("iML1515" => "E. coli K12",
"iJN678" => "Synechocystis\nsp. PCC6803",
"iJN746" => "P. putida",
"iRC1080" => "C. reinhardtii",
"iYO844" => "B. subtilis", 
"iMM904" => "S. cerevisiae", 
"iLJ478" => "T. maritima")

md = Dict{String, Float64}()
for modelname in models
    modelpath = joinpath("models", modelname) 
    model = CobraTools.read_model(modelpath)
    gibbs = CobraTools.map_gibbs_rxns(model.rxns) 
    nkeys = length([k for k in keys(gibbs) if !startswith(k, "EX_") && gibbs[k] != 0.0])
    nrxns = length(model.rxns)
    md[modelname[1:end-5]] = nkeys/nrxns
end

mdks = keys(md)
vs = [md[k] for k in mdks]
xnms = [k*"\n"*mnames[k] for k in mdks]
bar(vs)
plot!(xticks=(1:7,xnms), legend=false, ylabel="Fraction of reactions with ΔG⁰", xrotation=60)

