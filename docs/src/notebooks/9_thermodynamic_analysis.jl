
# # Thermodynamic analysis

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this notebook we will use max-min driving force analysis (MMDF) to find optimal concentrations
# for the metabolites in glycolysis to ensure that the smallest driving force across all the
# reactions in the model is as large as possible. For more information, see Flamholz, Avi, et al.
# "Glycolytic strategy as a tradeoff between energy yield and protein cost." Proceedings of
# the National Academy of Sciences 110.24 (2013): 10039-10044.

using COBREXA
using Tulip

# Create a model of glycolysis

#md # !!! tip "Directions of reactions"
#md #     Be careful when constructing models for MMDF, the reaction directions in the model
#md #     and thermodynamic data need to be consistent. You need to check each reaction to ensure
#md #     the ΔG'⁰ is in the direction of the reaction as saved in the model.

model = StandardModel("Glycolysis")

mets = [
    Metabolite("glc__D"),
    Metabolite("g6p"),
    Metabolite("f6p"),
    Metabolite("fdp"),
    Metabolite("dhap"),
    Metabolite("g3p"),
    Metabolite("13dpg"),
    Metabolite("3pg"),
    Metabolite("2pg"),
    Metabolite("pep"),
    Metabolite("pyr"),
    Metabolite("lac__D"),
    Metabolite("nadh"),
    Metabolite("nad"),
    Metabolite("h"),
    Metabolite("atp"),
    Metabolite("adp"),
    Metabolite("h2o"),
    Metabolite("pi"),
]

rxns = [
    Reaction(
        "HEX";
        metabolites = Dict(
            "atp" => -1.0,
            "glc__D" => -1.0,
            "g6p" => 1.0,
            "adp" => 1.0,
            "h" => 1.0,
        ),
    ),
    Reaction("PGI"; metabolites = Dict("g6p" => -1.0, "f6p" => 1.0)),
    Reaction(
        "PFK";
        metabolites = Dict(
            "f6p" => -1.0,
            "atp" => -1.0,
            "adp" => 1.0,
            "h" => 1.0,
            "fdp" => 1.0,
        ),
    ),
    Reaction("FBA"; metabolites = Dict("fdp" => -1.0, "dhap" => 1.0, "g3p" => 1.0)),
    Reaction("TPI"; metabolites = Dict("dhap" => -1.0, "g3p" => 1.0)),
    Reaction(
        "GAPD";
        metabolites = Dict(
            "g3p" => -1.0,
            "nad" => -1.0,
            "pi" => -1.0,
            "h" => 1.0,
            "nadh" => 1.0,
            "13dpg" => 1.0,
        ),
    ),
    Reaction(
        "PGK";
        metabolites = Dict("13dpg" => -1.0, "adp" => -1.0, "atp" => 1.0, "3pg" => 1.0),
    ),
    Reaction("PGM"; metabolites = Dict("3pg" => -1.0, "2pg" => 1)),
    Reaction("ENO"; metabolites = Dict("2pg" => -1.0, "h2o" => 1.0, "pep" => 1)),
    Reaction(
        "PYK";
        metabolites = Dict(
            "pep" => -1.0,
            "adp" => -1.0,
            "h" => -1.0,
            "atp" => 1.0,
            "pyr" => 1.0,
        ),
    ),
    Reaction(
        "LDH";
        metabolites = Dict(
            "pyr" => -1.0,
            "nadh" => -1.0,
            "h" => -1.0,
            "nad" => 1.0,
            "lac__D" => 1.0,
        ),
    ),
]

add_metabolites!(model, mets)
add_reactions!(model, rxns)

# Load some thermodynamic data, ΔG'⁰ from eQuilibrator

thermodynamic_data = Dict(
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

# Run max min driving force analysis with some reasonable constraints

df, dgs, concens = max_min_driving_force(
    model,
    Tulip.Optimizer,
    thermodynamic_data;
    proton_id = "h",
    water_id = "h2o",
    concentration_ratios = [("atp", "adp", 10.0), ("nadh", "nad", 0.1)],
    constant_concentrations = [("pi", 10e-3)],
    concentration_lb = 1e-6,
    concentration_ub = 10e-3,
)

# Plot results

relative_flux = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0] # scale fluxes
rids = [rxn.id for rxn in rxns]
dg_standard = [thermodynamic_data[rid] for rid in rids]
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
    (cumsum(dg_standard) .- first(dg_standard)) .* relative_flux;
    color = :red,
    label = "Standard",
)
lines!(
    ax,
    1:length(rids),
    (cumsum(dg_opt) .- first(dg_opt)) .* relative_flux;
    color = :blue,
    label = "Optimized",
)
ax.xticks = (1:length(rids), rids)
fig[1, 2] = Legend(fig, ax, "ΔG'", framevisible = false)
fig
