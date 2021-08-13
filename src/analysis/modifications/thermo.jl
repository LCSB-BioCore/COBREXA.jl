"""
    thermodynamic_constraints(
        gibbs_free_energies::Dict{String,Float64};
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

Add thermodynamic constraints to an optimization model. This function directly ensures that
the Gibbs energy dissipated by each reaction is negative, i.e. `ΔG * flux < 0`.
Consequently, quadratic constraints are added to the optimization model, causing the problem
to become nonlinear. This is unlike the traditional tFBA and llFBA that instead cast the
model into a mixed integer optimization problem. This modification adds `dgs` and `logcs` as
variables to the optimization model, where `dgs` are the actual Gibbs energy of each
reaction and `logcs` are the logged concentration of each metabolite. Only reactions and
metabolites involved with the supplied thermodynamic data have constraints placed on them,
all the others are free variables. This allows solvers to optimize these variables away, but
allows users to access the values of the solved problem in the order they are stored in the
model.

Since biological thermodynamics data are assumed to be passed into the function, protons and
water do not have adjustable concentration (constant pH and aqueous conditions are assumed).
Hence, they are always free variables. Additional constraints are supplied in
`concentration_ratios`, `constant_concentrations`, `concentration_lb`, and
`concentration_ub`. Their meaning is the same as in [`max_min_driving_force`](@ref).
Concentrations in Molar.

See also: [`thermodynamic_flux_balance_analysis_dict`](@ref) for a convenience wrapper of
this modification.
"""
thermodynamic_constraints(
    gibbs_free_energies::Dict{String,Float64};
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
) =
    (model, opt_model) -> begin

        # find reactions with thermodynamic data
        dg_rids = filter(x -> haskey(gibbs_free_energies, x), reactions(model))
        dg_ridxs = Int.(indexin(dg_rids, reactions(model)))

        #=
        remove protons, water from metabolites whose concentration can change, as well as all
        metabolites not involved in the reactions that have thermodynamic information
        =#
        dg_mids = unique(
            vcat(
                [
                    collect(keys(rxn.metabolites)) for
                    (rid, rxn) in model.reactions if rid in dg_rids
                ]...,
            ),
        )
        filter!(x -> !(x in [proton_ids; water_ids]), dg_mids)
        dg_midxs = Int.(indexin(dg_mids, metabolites(model)))

        # reduced stoichiometric matrix
        St = (stoichiometry(model)[dg_midxs, dg_ridxs])'

        RT = 298.15 * 8.314e-3 # kJ/mol

        dg0s = [gibbs_free_energies[rid] for rid in dg_rids] # NB: in order of dg_ridxs

        # Add variables
        @variables opt_model begin
            dgs[1:n_reactions(model)] # all reactions get dG data, but only the reactions with data are constrained
            logcs[1:n_metabolites(model)] # all metabolites become variables - solver will remove unnecessary ones
        end

        @constraints opt_model begin
            dgs[dg_ridxs] .== dg0s .+ RT .* St * logcs[dg_midxs]
            dgs[dg_ridxs] .* opt_model[:x][dg_ridxs] .<= 0.0 # thermodynamically consistent
        end

        log_lb = log(concentration_lb)
        log_ub = log(concentration_ub)
        constant_concentration_dict =
            Dict(met => val for (met, val) in constant_concentrations)
        # add bounds for metabolites, metabolites not involved in thermo as free
        for mid in dg_mids
            i = first(indexin([mid], metabolites(model)))
            if haskey(constant_concentration_dict, mid)
                @constraint(opt_model, logcs[i] == log(constant_concentration_dict[mid]))
            else
                @constraint(opt_model, log_lb <= logcs[i] <= log_ub)
            end
        end

        # add metabolite ratio constraints
        for (met1, met2, val) in concentration_ratios
            i, j = Int.(indexin([met1, met2], metabolites(model)))
            @constraint(opt_model, logcs[i] == log(val) + logcs[j])
        end
    end
