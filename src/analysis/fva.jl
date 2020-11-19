"""
Flux variability analysis (FVA)

FVA solves the pair of optimization problems for each flux xᵢ
```
min/max xᵢ
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
     cᵀx ≥ γ Z₀
```
where Z₀:= cᵀx₀ is the objective value of an optimal solution to the associated
FBA problem
"""
function fluxVariabilityAnalysis(model::LinearModel, optimizer)
   (m, n) = size(model.S)
   return fluxVariabilityAnalysis(model, collect(1:n), optimizer)
end

function fluxVariabilityAnalysis(model::LinearModel, reactions::Array{Int64, 1}, optimizer)
   (maximum(reactions) > length(model.rxns)) && error("Index exceeds number of reactions.")
   γ = 1.
   fluxes = zeros(length(reactions), 2)

   (optimization_model, x₀) = fluxBalanceAnalysis(model::LinearModel, optimizer)
   Z₀ = JuMP.objective_value(optimization_model)
   x = all_variables(optimization_model)
   @constraint(optimization_model, model.c' * x ≥ γ * Z₀)

   for i in eachindex(reactions)
      sense = MOI.MIN_SENSE
      @objective(optimization_model, sense, x[reactions[i]])
      JuMP.optimize!(optimization_model)
      fluxes[i, 1] = JuMP.objective_value(optimization_model)

      sense = MOI.MAX_SENSE
      JuMP.set_objective_sense(optimization_model, sense)
      JuMP.optimize!(optimization_model)
      fluxes[i, 2] = JuMP.objective_value(optimization_model)
   end
   return fluxes
end

"""
(Local) multi-process version of FVA that assigns the reactions to several
processes running in parallel
`workers` should be a list of process ids returned by `createParPool`
"""
function parFVA(model::LinearModel, reactions::Array{Int64, 1}, solverName::Symbol, workersToUse::Array{Int64,1})
   nReacs = length(reactions)
   nWorkers = length(workersToUse)
   if nReacs < nWorkers
      @info "Number of workers exceeds number of reactions. 1 worker per reaction will be used."
      workersToUse = workersToUse[1:nReacs]
      nWorkers = nReacs
   end
   remRefs = Array{Future}(undef, nWorkers)
   fluxes = zeros(nReacs, 2)
   alloReacs = allocateReacs(reactions, nWorkers)

   @sync for (round, pid) in enumerate(workersToUse)
      @async remRefs[round] = @spawnat pid begin
         fluxVariabilityAnalysis(model, alloReacs[round], eval(solverName).Optimizer)
      end
   end
   @sync for (round, pid) in enumerate(workersToUse)
      fluxes[alloReacs[round], :] = fetch(remRefs[round])
   end

   return fluxes
end

"""
Auxiliary function for `parFVA` to divide the list of reactions evenly (for now)
"""
function allocateReacs(reactions::Array{Int64, 1}, nWorkers::Int)
   nReacs = length(reactions)
   steps = floor(Int, nReacs/nWorkers) * ones(Int, nWorkers)
   steps[1:nReacs%nWorkers] = steps[1:nReacs%nWorkers] .+ 1
   steps = cumsum(steps)
   allocatedReacs = Vector(undef, nWorkers)

   startI = 1
   for (i, endI) in enumerate(steps)
      allocatedReacs[i] = reactions[startI:endI]
      startI = endI + 1
   end

   return allocatedReacs
end
