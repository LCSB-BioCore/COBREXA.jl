"""
$(TYPEDSIGNATURES)

Change the abundances of the community model within the solver. Only
[`CommunityModel`](@ref) and [`EqualGrowthCommunityModel`](@ref) are supported
at this time.
"""
modify_abundances(new_abundances::Vector{Float64}) =
    (model, opt_model) -> begin
        #=
        Only support these community models because the location of the
        environmental balances are known.
        =#
        model isa CommunityModel ||
            model isa EqualGrowthCommunityModel ||
            throw(
                ArgumentError(
                    "Only CommunityModel and EqualGrowthCommunityModel are supported at this time.",
                ),
            )

        check_abundances(new_abundances) # TODO consider removing: too pedantic

        env_rows =
            model isa CommunityModel ?
            environment_exchange_stoichiometry(model, new_abundances) :
            environment_exchange_stoichiometry(model.inner, new_abundances)
        env_link = spdiagm(sum(env_rows, dims = 2)[:])

        n_vars = n_variables(model)
        n_env_vars = length(model.environmental_links)
        n_cons = length(opt_model[:mb])
        n_objs = model isa CommunityModel ? 0 : length(model.inner.members)

        row_offset =
            model isa CommunityModel ? n_cons - n_env_vars : n_cons - n_env_vars - n_objs

        # fix abundance coefficients of species exchanges
        for (i, j, v) in zip(findnz(env_rows)...)
            ii = i + row_offset
            set_normalized_coefficient(opt_model[:mb][ii], opt_model[:x][j], v)
        end

        column_offset =
            model isa CommunityModel ? n_vars - n_env_vars : n_vars - n_env_vars - 1

        # fix total abundance to link exchange
        for (i, j, v) in zip(findnz(env_link)...)
            jj = j + column_offset
            ii = i + row_offset
            set_normalized_coefficient(opt_model[:mb][ii], opt_model[:x][jj], -v)
        end
    end
