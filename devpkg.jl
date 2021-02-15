using CobraTools
using Measurements
using LinearAlgebra
using JLD
using JuMP
using Gurobi
using SparseArrays


# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.read_model(modelpath)

biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
fbasol = CobraTools.fba(model, biomass_rxn)
pfbasol = CobraTools.pfba(model, biomass_rxn)
ad = CobraTools.atom_exchange(pfbasol)


# ecoli_kJmol = -71.36 ± 8.55 # formation of biomass kJ/mol
# gibbs = CobraTools.map_gibbs_rxns(model.rxns, dgtype="prime") # slow 
# CobraTools.saveGibbs(joinpath("data", "dgprime.jld"), gibbs)
# gibbs = CobraTools.load_Gibbs(joinpath("data", "dgzeros.jld")) 
# CobraTools.map_gibbs_external(pfbasol, gibbs)
# CobraTools.map_gibbs_internal(pfbasol, gibbs)

cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)
set_optimizer(cbmodel, Gurobi.Optimizer)
set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]


# Fix glucose use 1.0 then normalization is easy. NB - if not 1 then change normalization!!
CobraTools.set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)

# Aerobic
# CobraTools.set_bound(o2_index, ubs, lbs; ub=1000.0, lb=-1000.0)
# Anaerobic
CobraTools.set_bound(o2_index, ubs, lbs; ub=1000.0, lb=0.0)

# No free ATP generation
CobraTools.set_bound(atpm_index, ubs, lbs; ub=1000.0, lb=0.0)

# FBA
@objective(cbmodel, Max, v[biomass_index])
optimize!(cbmodel) 
termination_status(cbmodel) != MOI.OPTIMAL && @warn "Optimization issue..."
μ_max = round(objective_value(cbmodel), digits=6)

CobraTools.set_bound(biomass_index, ubs, lbs; ub=μ_max, lb=μ_max)

wpoints, var_rxn_inds = CobraTools.get_warmup_points(cbmodel, v, ubs, lbs)
# samples = CobraTools.achr(10, warmup_points, var_rxn_inds, ubs, lbs)

###
N = 10
# Minimum allowed distance to the closest constraint
maxMinTol = 1e-9
# Ignore directions where u is really small
uTol = 1e-9
# Project out of directions that are too close to the boundary
dTol = 1e-14


N = nullspace(Array(S))
P = N * inv(transpose(N)*N) * transpose(N) # project v to ̂v that is satisfies the mass balances  P*v => ̂v    

nwpts = size(wpoints, 2)
center_point = mean(wpoints, dims=2)

init_index = rand(1:nwpts)
samples = zeros(size(wpoints, 1), N)
samples[:, 1] .= wpoints[:, init_index]

lbs = zeros(length(lbconref))
for i in eachindex(lbs)
    lbval = normalized_rhs(lbconref[i])
    if lbval > 0
        lbs[i] = -abs(lbval)
    else
        lbs[i] = abs(lbval)
    end
end
ubs = [normalized_rhs(ubconref[i]) for i in eachindex(ubconref)]

u = zeros(size(wpoints, 1)) # direction vector
λ = 0.0 # step size
dist_ubs = zeros(length(ubs))
dist_lbs = zeros(length(lbs))
# for n in 2:N
n=2    
rand_index1 = rand(1:nwpts)
    rand_index2 = rand(1:n-1)
    u .= rand() <= nwpts/(nwpts+n-1) ? wpoints[:, rand_index1] : samples[:, rand_index2]
    u .= u./norm(u)

    dist_lbs .= samples[:, n] - lbs
    dist_ubs .= ubs - samples[:, n]

    valid_dirs = (dist_ubs .> dTol) .== (dist_lbs .> dTol)
    
    posDirn = u[valid_dirs] .> uTol
    negDirn = u[valid_dirs] .< -uTol

    maxStepTemp = dist_ubs[valid_dirs]./u[valid_dirs]
    minStepTemp = -dist_lbs[valid_dirs]./u[valid_dirs]
    maxStepVec = [maxStepTemp[posDirn]; minStepTemp[negDirn]]
    minStepVec = [minStepTemp[posDirn]; maxStepTemp[negDirn]]
    
    maxStep = minimum(maxStepVec)
    minStep = maximum(minStepVec)

    if (abs(minStep) < maxMinTol && abs(maxStep) < maxMinTol) || (minStep > maxStep)
        println("Warning $(minStep) $(maxStep)")
        # continue
    end

    stepDist = rand()*(maxStep-minStep)+minStep

    # Advance to the next point
    curPoint = prevPoint + stepDist*u;

# end















# # Yeast GEM
# modelpath = joinpath("models", "iMM904.json") 
# model = CobraTools.readmodel(modelpath)

# biomass_rxn = findfirst(model.rxns, "BIOMASS_SC5_notrace")
# fbasol = CobraTools.fba(model, biomass_rxn)
# pfbasol = CobraTools.pfba(model, biomass_rxn)

# ad = CobraTools.atom_exchange(pfbasol)
