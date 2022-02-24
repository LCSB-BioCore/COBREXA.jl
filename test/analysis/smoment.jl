using COBREXA,
    Tulip,
    JSON,
    JuMP,
    Statistics,
    TerminalPager,
    LinearAlgebra,
    ForwardDiff,
    CairoMakie,
    SparseArrays


#! remove isozymes with lower effectivity (kcat/total_mass), only one enzyme per reaction
remove_slow_isozymes!(model, reaction_kcats, protein_stoichiometry, protein_masses)

#: SMOMENT
total_protein_mass = 200.0 # mg/gDW
model.reactions["EX_glc__D_e"].lb = -1000.0 #! unconstrain otherwise bound will be hit
obj_id = "BIOMASS_Ecoli_core_w_GAM§FOR"

c, E, d, M, h, reaction_map, metabolite_map = smoment(
    model,
    protein_stoichiometry,
    protein_masses,
    reaction_kcats;
    total_protein_mass,
);

m = Model(Gurobi.Optimizer);
@variable(m, x[1:size(E, 2)]);
bid = reaction_map[obj_id]
@objective(m, Max, x[bid]);
@constraint(m, E * x .== d);
@constraint(m, M * x .<= h);
optimize!(m)

reaction_fluxes_pre = map_ids_to_sols(reaction_map, value.(x));
plot_flux_summary(reaction_fluxes_pre)