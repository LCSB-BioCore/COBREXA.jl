"""
    thermodynamic_flux_balance_analysis_dict(model, optimizer; kwargs...)

A shortcut for wrapping the output of [`flux_balance_analysis`](@ref) with thermodyamic
data, `kwargs` must contain `modifications` where [`thermodynamic_constraints`](@ref) is
included for this convenience wrapper to make sense. Returns a named tuple mapping `fluxes`,
`concentrations` and `gibbs_reaction_energies` to reactions and metabolites in `model`.

# Example
```
sol = thermodynamic_flux_balance_analysis_dict(
    model,
    optimizer;
    modifications = [
        thermodynamic_constraints(
            gibbs_free_energies;
            proton_ids = ["h_c", "h_e"],
            water_ids = ["h2o_c", "h2o_e"],
            concentration_ratios = [
                ("atp_c", "adp_c", 10.0),
                ("adp_c", "amp_c", 1.0),
                ("nadph_c", "nadp_c", 10.0),
                ("nadh_c", "nad_c", 0.1),
            ],
            constant_concentrations = [
                ("coa_c", 1e-3),
                ("co2_c", 10e-6),
                ("pi_c", 10e-3),
                ("ppi_c", 1e-3),
            ],
            concentration_lb = 1e-6,
            concentration_ub = 10e-3,
        )
    ],
)
```
"""
function thermodynamic_flux_balance_analysis_dict(args...; kwargs...)

    sol = flux_balance_analysis(args...; kwargs...)

    isnothing(sol) && return nothing

    fluxes = Dict(zip(reactions(args[1]), value.(sol[:x])))
    concens = Dict(zip(metabolites(args[1]), exp.(value.(sol[:logcs]))))
    dgs = Dict(zip(reactions(args[1]), value.(sol[:dgs])))
    return (fluxes = fluxes, concentrations = concens, gibbs_reaction_energies = dgs)
end
