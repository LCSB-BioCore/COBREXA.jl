using COBREXA, Gurobi, CSV, JSONFBCModels, JSON
import AbstractFBCModels as A
import ConstraintTrees as C
using CairoMakie

# get model
download_model(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    "iML1515.json",
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

model = load_model("iML1515.json")

# get reaction kcats from Heckmann 2020 Nature Comms
rid_kcats = Dict(
    x.react_id => x.kappmax_KO_ALE_davidi_per_pp_per_s_ensemble_model for
    x in CSV.File(joinpath("docs", "src", "showcase", "Dataset_S1C_turnonver_n.csv"))
) # 1/s

reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())
        kcat_f = get(rid_kcats, rid, 30.0)
        kcat_r = get(rid_kcats, rid * "_b", kcat_f)
        d["isozyme_"*string(i)] = Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = kcat_f * 3.6, # forward reaction turnover number units = k/h
            kcat_reverse = kcat_r * 3.6, # reverse reaction turnover number units = k/h
        )
    end
end

# from Uniprot
gene_product_molar_masses = Dict(
    k => v for
    (k, v) in JSON.parsefile(joinpath("docs", "src", "showcase", "ecoli_gpmms.tsv"))
)
gene_product_molar_masses["b1692"] = 54.58 # missing
gene_product_molar_masses["s0001"] = 54.58 # missing

# wt simulation
wt = flux_balance_constraints(model; interface = :identifier_prefixes)
wt.fluxes.EX_glc__D_e.bound = C.Between(-1000, 1000)
ec_wt = add_enzyme_constraints(
    wt;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 400.0, # mg/gDW
)

wt_sol = optimized_constraints(
    ec_wt;
    optimizer = Gurobi.Optimizer,
    objective = ec_wt.objective.value,
    settings = [silence,]
)


# KO simulations
kos = Dict(
    "argA" => (:ACGS, :arg),
    "serA" => (:PGCD, :ser),
    "ilvA" => (:THRD_L, :ile),
    "metA" => (:HSST, :met),
    "hisB" => (:HISTP, :his),
    "proA" => (:G5SD, :pro),
    "cysE" => (:SERAT, :cys),
    "thrC" => (:THRS, :thr),
    "glyA" => (:GHMT2r, :gly),
    "leuB" => (:IPMD, :leu),
    "trpC" => (:IGPS, :trp),
    "pheA" => (:PPNDH, :phe),
    "tyrA" => (:PPND, :tyr),
    "lysA" => (:DAPDC, :lys),
)


open_rxns = [
    :EX_glu__L_e,
]

measured_ab1 = 0.59
ko1 = "cysE"
ko2 = "tyrA"
(ko_rxn1, aa1) = kos[ko1]
(ko_rxn2, aa2) = kos[ko2]

aa1_id = aa1 == :gly ? Symbol("EX_", aa1, "_e") : Symbol("EX_", aa1, "__L_e")
aa2_id = aa2 == :gly ? Symbol("EX_", aa2, "_e") : Symbol("EX_", aa2, "__L_e")

aa1_flux = abs(wt_sol.fluxes[ko_rxn1])*5
aa2_flux = abs(wt_sol.fluxes[ko_rxn2])*5

# create a isoleucine knockout 
mut1 = flux_balance_constraints(model; interface = :identifier_prefixes)
mut1.fluxes.EX_glc__D_e.bound = C.Between(-1000, 1000)
mut1.fluxes[ko_rxn1].bound = C.EqualTo(0.0) # knockout

for (_, aa) in values(kos)
    aa_id = aa == :gly ? Symbol("EX_", aa, "_e") : Symbol("EX_", aa, "__L_e")
    mut1.fluxes[aa_id].bound = C.Between(-1000, 1000)    
end

# mut1.fluxes[aa1_id].bound = C.Between(-aa1_flux, aa1_flux)
# mut1.fluxes[aa2_id].bound = C.Between(-aa2_flux, aa2_flux)

ec_mut1 = add_enzyme_constraints(
    mut1;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 400.0, # mg/gDW
)

# create a methionine knockout
mut2 = flux_balance_constraints(model; interface = :identifier_prefixes)
mut2.fluxes.EX_glc__D_e.bound = C.Between(-1000, 1000)
mut2.fluxes[ko_rxn2].bound = C.EqualTo(0.0) # knockout

for (_, aa) in values(kos)
    aa_id = aa == :gly ? Symbol("EX_", aa, "_e") : Symbol("EX_", aa, "__L_e")
    mut2.fluxes[aa_id].bound = C.Between(-1000, 1000)    
end

# mut2.fluxes[aa1_id].bound = C.Between(-aa1_flux, aa1_flux)
# mut2.fluxes[aa2_id].bound = C.Between(-aa2_flux, aa2_flux)

ec_mut2 = add_enzyme_constraints(
    mut2;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 400.0, # mg/gDW
)

# create bound function to constrain interface
boundf(id) = begin
    ex_id = first(id)
    if ex_id == aa1_id || ex_id == aa2_id
        C.EqualTo(0.0)
    else
        mut1.interface.exchanges[ex_id].bound # have same interface, so easy
    end
end

growthrates = Float64[]
abundances = collect(range(0, 1, length=100))
for a in abundances
    x = interface_constraints(
        :mut1 => (ec_mut1, ec_mut1.interface.exchanges, a),
        :mut2 => (ec_mut2, ec_mut2.interface.exchanges, 1-a);
        bound = boundf
    )
    
    x *=
        :equalgrowth^C.Constraint(
            C.value(x.mut1.fluxes.BIOMASS_Ec_iML1515_core_75p37M) -
            C.value(x.mut2.fluxes.BIOMASS_Ec_iML1515_core_75p37M),
            C.EqualTo(0.0),
        )
    
    sol = optimized_constraints(
        x;
        optimizer = Gurobi.Optimizer,
        objective = x.mut1.objective.value,
        settings = [silence,]
    )

    push!(growthrates, sol.mut1.objective)
end


fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="Abundance",
    ylabel="Growth rate",
    title="$(ko1) + $(ko2)"
)
lines!(ax, abundances, growthrates)
vspan!(ax, measured_ab1-0.05, measured_ab1+0.05)
fig

