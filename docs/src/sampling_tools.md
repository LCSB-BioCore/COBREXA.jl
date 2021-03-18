# Sampling Tools
It is well known that the FBA does not yield a unique solution, i.e. many flux distributions are capable of satisfying the system constraints as well as optimizing the imposed objective function. 
Let the feasible space be defined by ``\mathcal{P} = \left\{ v : Sv = 0 \cap v_{\text{min}} \leq v \leq v_{\text{max}} \right\}``. 
Sampling methods have been developed to uniformly sample from this feasible solution space. 
`COBREXA.jl` implements both `hit_and_run` and `achr` to sample from ``\mathcal{P}``.
```@docs
hit_and_run
achr
``` 
```@example sample
using COBREXA
using JuMP
using Tulip

download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")
model = read_model("e_coli_core.json")
rm("e_coli_core.json") # hide
optimizer = Tulip.Optimizer
biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
cons = Dict("EX_glc__D_e" => (-12.0, -12.0))
sol = fba(model, biomass, optimizer, constraints=cons) # classic flux balance analysis
cons["BIOMASS_Ecoli_core_w_GAM"] = (sol["BIOMASS_Ecoli_core_w_GAM"], sol["BIOMASS_Ecoli_core_w_GAM"]*0.99)

samples = hit_and_run(100_000, model, optimizer; keepevery=10, samplesize=5000, constraints=cons)
```
