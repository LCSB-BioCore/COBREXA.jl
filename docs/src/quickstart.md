
# Quick start guide

You can install COBREXA from Julia repositories. Start `julia`, **press `]`** to
switch to the Packaging environment, and type:
```
add COBREXA
```

You also need to install your favorite solver supported by `JuMP.jl` (such as
Gurobi, Mosek, CPLEX, GLPK, Clarabel, etc., see a [list
here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)).  For
example, you can install `Tulip.jl` solver by typing:
```
add Tulip
```

Alternatively, you may use [prebuilt Docker and Apptainer images](#prebuilt-images).

If you are running COBREXA.jl for the first time, it is very likely that upon
installing and importing the packages, your Julia installation will need to
precompile their source code from the scratch. In fresh installations, the
precompilation process should take less than 5 minutes.

When the packages are installed, switch back to the "normal" julia shell by
pressing Backspace (the prompt should change color back to green). After that,
you can download [a SBML model from the
internet](http://bigg.ucsd.edu/models/e_coli_core) and perform a
flux balance analysis as follows:

```julia
using COBREXA   # loads the package
using Tulip     # loads the optimization solver

# download the model
download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

# open the SBML file and load the contents
model = load_model("e_coli_core.xml")

# run a FBA
fluxes = flux_balance_analysis_dict(model, Tulip.Optimizer)
```

The variable `fluxes` will now contain a dictionary of the computed optimal
flux of each reaction in the model:
```
Dict{String,Float64} with 95 entries:
  "R_EX_fum_e"    => 0.0
  "R_ACONTb"      => 6.00725
  "R_TPI"         => 7.47738
  "R_SUCOAS"      => -5.06438
  "R_GLNS"        => 0.223462
  "R_EX_pi_e"     => -3.2149
  "R_PPC"         => 2.50431
  "R_O2t"         => 21.7995
  "R_G6PDH2r"     => 4.95999
  "R_TALA"        => 1.49698
  ⋮               => ⋮
```

#### Model variant processing

The main feature of COBREXA.jl is the ability to easily specify and process
many analyses in parallel. To demonstrate, let's see how the organism would
perform if some reactions were disabled independently:

```julia
# convert to a model type that is efficient to modify
m = convert(ObjectModel, model)

# find the model objective value if oxygen or carbon dioxide transports are disabled
screen(m, # the base model
    variants=[ # this specifies how to generate the desired model variants
        [], # one with no modifications, i.e. the base case
        [with_changed_bound("R_O2t", lower_bound =0.0, upper_bound =0.0)], # disable oxygen
        [with_changed_bound("R_CO2t", lower_bound =0.0, upper_bound =0.0)], # disable CO2
        [with_changed_bound("R_O2t", lower_bound =0.0, upper_bound =0.0),
	        with_changed_bound("R_CO2t", lower_bound =0.0, upper_bound =0.0)], # disable both
    ],
    # this specifies what to do with the model variants (received as the argument `x`)
    analysis = x ->
        flux_balance_analysis_dict(x, Tulip.Optimizer)["R_BIOMASS_Ecoli_core_w_GAM"],
)
```
You should receive a result showing that missing oxygen transport makes the
biomass production much harder:
```julia
4-element Vector{Float64}:
 0.8739215022674809
 0.21166294973372796
 0.46166961413944896
 0.21114065173865457
```

Most importantly, such analyses can be easily specified by automatically
generating long lists of modifications to be applied to the model, and
parallelized.

Knocking out each reaction in the model is efficiently accomplished:

```julia
# load the task distribution package, add several worker nodes, and load
# COBREXA and the solver on the nodes
using Distributed
addprocs(4)
@everywhere using COBREXA, Tulip

# get a list of the workers
worker_list = workers()

# run the processing in parallel for many model variants
res = screen(m,
    variants=[
	# create one variant for each reaction in the model, with that reaction knocked out
        [with_changed_bound(reaction_id, lower_bound =0.0, upper_bound =0.0)]
	for reaction_id in reactions(m)
    ],
    analysis = model -> begin
	# we need to check if the optimizer even found a feasible solution,
	# which may not be the case if we knock out important reactions
    	sol = flux_balance_analysis_dict(model, Tulip.Optimizer)
	isnothing(sol) ? nothing : sol["R_BIOMASS_Ecoli_core_w_GAM"]
    end,
    # run the screening in parallel on all workers in the list
    workers = worker_list,
)
```

In result, you should get a long list of the biomass production for each
reaction knockout. Let's decorate it with reaction names:
```julia
Dict(reactions(m) .=> res)
```
...which should output an easily accessible dictionary with all the objective
values named, giving a quick overview of which reactions are critical for the
model organism to create biomass:
```julia
Dict{String, Union{Nothing, Float64}} with 95 entries:
  "R_ACALD"       => 0.873922
  "R_PTAr"        => 0.873922
  "R_ALCD2x"      => 0.873922
  "R_PDH"         => 0.796696
  "R_PYK"         => 0.864926
  "R_CO2t"        => 0.46167
  "R_EX_nh4_e"    => 1.44677e-15
  "R_MALt2_2"     => 0.873922
  "R_CS"          => 2.44779e-14
  "R_PGM"         => 1.04221e-15
  "R_TKT1"        => 0.864759
  ⋮             => ⋮
```
