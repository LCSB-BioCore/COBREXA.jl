
# # Maximum-minimum driving force analysis

#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this notebook we use the max-min driving force analysis (MMDFA) to find
# optimal concentrations for the metabolites in glycolysis to ensure that the
# smallest driving force across all the reactions in the model is as large as
# possible. For more information, see Flamholz, Avi, et al.  "Glycolytic
# strategy as a tradeoff between energy yield and protein cost.", Proceedings
# of the National Academy of Sciences 110.24 (2013): 10039-10044.

using COBREXA, Tulip

# Let's first make a model of glycolysis fermentation.

mets = [
    "13dpg",
    "2pg",
    "3pg",
    "adp",
    "atp",
    "dhap",
    "f6p",
    "fdp",
    "g3p",
    "g6p",
    "glc__D",
    "h",
    "h2o",
    "lac__D",
    "nad",
    "nadh",
    "pep",
    "pi",
    "pyr",
]

rxns = Dict(
    "ENO" => Dict("2pg" => -1.0, "h2o" => 1.0, "pep" => 1.0),
    "FBA" => Dict("fdp" => -1.0, "dhap" => 1.0, "g3p" => 1.0),
    "GAPD" => Dict(
        "g3p" => -1.0,
        "nad" => -1.0,
        "pi" => -1.0,
        "h" => 1.0,
        "nadh" => 1.0,
        "13dpg" => 1.0,
    ),
    "HEX" =>
        Dict("atp" => -1.0, "glc__D" => -1.0, "g6p" => 1.0, "adp" => 1.0, "h" => 1.0),
    "LDH" =>
        Dict("pyr" => -1.0, "nadh" => -1.0, "h" => -1.0, "nad" => 1.0, "lac__D" => 1.0),
    "PFK" => Dict("f6p" => -1.0, "atp" => -1.0, "adp" => 1.0, "h" => 1.0, "fdp" => 1.0),
    "PGI" => Dict("g6p" => -1.0, "f6p" => 1.0),
    "PGK" => Dict("13dpg" => -1.0, "adp" => -1.0, "atp" => 1.0, "3pg" => 1.0),
    "PGM" => Dict("3pg" => -1.0, "2pg" => 1.0),
    "PYK" =>
        Dict("pep" => -1.0, "adp" => -1.0, "h" => -1.0, "atp" => 1.0, "pyr" => 1.0),
    "TPI" => Dict("dhap" => -1.0, "g3p" => 1.0),
)

model = StandardModel("Glycolysis")

add_metabolites!(model, Metabolite.(mets))
add_reactions!(model, collect(Reaction(rid; metabolites = mets) for (rid, mets) in rxns))

model

# We need some thermodynamic data. You can get Gibbs free energies (ΔG⁰) e.g.
# from [eQuilibrator](https://equilibrator.weizmann.ac.il/), possibly using the
# [Julia wrapper](https://github.com/stelmo/eQuilibrator.jl) that allows you to
# automate this step. Here, we make a dictionary that maps the reaction IDs to
# calculated Gibbs free energies of reactions.

gibbs_free_energies = Dict( # ΔG⁰ in kJ/mol
    "TPI" => 5.57535,
    "PGK" => -19.32,
    "PFK" => -14.5988,
    "ENO" => -3.81089,
    "PYK" => -27.5833,
    "LDH" => -23.6803,
    "FBA" => 22.3932,
    "PGI" => 2.6617,
    "GAPD" => 4.60271,
    "PGM" => -4.52041,
    "HEX" => -17.90,
)

# Run max min driving force analysis with some reasonable constraints. Protons
# and water are removed from the concentration calculation of the optimization
# problem, thus we specify their IDs explicitly.  The reason for this is that
# the standard Gibbs free energy change of biochemical reactions take place at
# constant pH, so proton concentration should not change to make the analysis
# behave reasonably; likewise we just assume that reactions occur in relatively
# stable aqueous environments, hence water excluded too.

df, dgs, concens = max_min_driving_force(
    model,
    gibbs_free_energies,
    Tulip.Optimizer;
    ignore_metabolites = ["h", "h2o"],
    concentration_ratios = Dict(("atp", "adp") => 10.0, ("nadh", "nad") => 0.1),
    constant_concentrations = Dict("pi" => 1e-2), # constant phosphate concentration set to 10 mM
    concentration_lb = 1e-6, # minimum 1 μM for all metabolites
    concentration_ub = 1e-2, # maximum 10 mM or all metabolites
)

# Plot the results to show how the concentrations can be used to ensure that
# each reach proceeds "down hill" (ΔᵣG < 0) and that the driving force is as
# large as possible across all the reactions in the model. Compare this to the
# driving forces at standard conditions.

# We additionally scale the fluxes according to their stoichiometry in the
# pathway. From the output, it is clear that that metabolite concentrations
# play a large role in ensuring the thermodynamic consistency of in vivo enzyme
# reactions.

relative_flux = Dict(
    "HEX" => 1.0,
    "PGI" => 1.0,
    "PFK" => 1.0,
    "FBA" => 1.0,
    "TPI" => 1.0,
    "GAPD" => 2.0,
    "PGK" => 2.0,
    "PGM" => 2.0,
    "ENO" => 2.0,
    "PYK" => 2.0,
    "LDH" => 2.0,
)

rids = ["HEX", "PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK", "LDH"] # in order of pathway
rid_rf = [relative_flux[rid] for rid in rids]
dg_standard = [gibbs_free_energies[rid] for rid in rids]
dg_opt = [dgs[rid] for rid in rids]

using CairoMakie

fig = Figure();
ax = Axis(
    fig[1, 1],
    xticklabelrotation = -pi / 2,
    xlabel = "Reaction",
    ylabel = "Cumulative ΔG [kJ/mol]",
);

lines!(
    ax,
    1:length(rids),
    (cumsum(dg_standard) .- first(dg_standard)) .* rid_rf;
    color = :red,
    label = "Standard",
)
lines!(
    ax,
    1:length(rids),
    (cumsum(dg_opt) .- first(dg_opt)) .* rid_rf;
    color = :blue,
    label = "Optimized",
)
ax.xticks = (1:length(rids), rids)
fig[1, 2] = Legend(fig, ax, "ΔG'", framevisible = false)
fig

#md # !!! tip "Directions of reactions"
#md #     Be careful when constructing models for MMDFA, the reaction directions in the model
#md #     and thermodynamic data need to be consistent with the overall flux
#md #     direction implied by the model. For example, in BiGG, `LDH_D` is written
#md #     `lac__D + nad ⟷ h + nadh + pyr` and the associated ΔrG'⁰ is 23.6803 kJ/mol.
#md #     For MMDFA no flux is calculated, so you need to write the reaction
#md #     in the direction of flux, i.e. `h + nadh + pyr ⟶ lac__D + nad` with ΔrG'⁰ as
#md #     -23.6803 kJ/mol.
