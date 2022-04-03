"""
    add_moment_constraints(
        ksas::Dict{String,Float64},
        protein_mass_fraction::Float64
    )

A modification that adds enzyme capacity constraints to the problem using a _modified_
version of the MOMENT algorithm. Requires specific activities, `ksas` [mmol product/g
enzyme/h], for each reaction. Proteins are identified by their associated gene IDs. Adds a
variable vector `y` to the problem corresponding to the protein concentration [g enzyme/gDW
cell] of each gene product in the order of `genes(model)`. The total protein concentration
[g protein/gDW cell] is constrained to be less than or equal to the `protein_mass_fraction`.
Reaction flux constraints are changed to the MOMENT constraints (see below) for all
reactions that have a gene reaction rule, otherwise the flux bounds are left unaltered.

See Adadi, Roi, et al. "Prediction of microbial growth rate versus biomass yield by a
metabolic network with kinetic parameters." PLoS computational biology (2012) for more
details of the original algorithm.

Here, a streamlined version of the algorithm is implemented to ensure that the correct units
are used. Specifically, this implementation uses specific activities instead of `kcats`.
Thus, for a reaction that can only proceed forward and is catalyzed by protein `a`, the flux
`x[i]` is bounded by `x[i] <= ksas[i] * y[a]`. If isozymes `a` or `b` catalyse the
reaction, then `x[i] <= ksas[i] * (y[a] + y[b])`. If a reaction is catalyzed by subunits `a`
and `b` then `x[i] <= ksas[i] * min(y[a], y[b])`. These rules are applied recursively in the
model like in the original algorithm. The enzyme capacity constraint is then implemented by
`sum(y) ≤ protein_mass_fraction`. The major benefit of using `ksas` instead of `kcats` is
that active site number and unit issues are prevented.

# Example
```
flux_balance_analysis(
    ...,
    modifications = [ add_moment_constraints(my_kcats, 0.6) ],
)
"""
add_moment_constraints(kcats::Dict{String,Float64}, protein_mass_fraction::Float64) =
    (model, opt_model) -> begin
        @warn(
            "DEPRECATION WARNING: 'add_moment_constraints' will be removed in future versions of COBREXA.jl in favor of a GECKO-based formulation"
        )

        lbs, ubs = get_optmodel_bounds(opt_model) # to assign directions
        # get grrs and ignore empty blocks: TODO: fix importing to avoid this ugly conditional see #462
        grrs = Dict(
            rid => reaction_gene_association(model, rid) for
            rid in reactions(model) if !(
                reaction_gene_association(model, rid) == [[""]] ||
                isnothing(reaction_gene_association(model, rid))
            )
        )

        # add protein variables
        y = @variable(opt_model, y[1:n_genes(model)] >= 0)

        # add variables to deal with enzyme subunits
        num_temp = sum(length(values(grr)) for grr in values(grrs)) # number "AND" blocks in grrs
        t = @variable(opt_model, [1:num_temp]) # anonymous variable
        @constraint(opt_model, t .>= 0)

        #=
        Note, not all of t needs to be created, only those with OR GRR rules, however this
        adds a lot of complexity to the code - re-implement this method if efficiency becomes
        an issue.
        =#

        # add capacity constraint
        @constraint(opt_model, sum(y) <= protein_mass_fraction)

        k = 1 # counter
        kstart = 1
        kend = 1
        x = opt_model[:x]
        for (ridx, rid) in enumerate(reactions(model))

            isnothing(get(grrs, rid, nothing)) && continue # only apply MOMENT constraints to reactions with a valid GRR
            grrs[rid] == [[""]] && continue # TODO: remove once #462 is implemented

            # delete original constraints
            delete(opt_model, opt_model[:lbs][ridx])
            delete(opt_model, opt_model[:ubs][ridx])

            #=
            For multi-subunit enzymes, the flux is bounded by the minimum concentration of
            any subunit in the enzyme. If multiple isozymes with subunits exist, then the
            sum of these minima bound the enzyme flux. E.g., suppose that you have a GRR
            [[y, z], [w, u, v]], then x <= sum(min(y, z), min(w, u, v)).

            It is possible to reformulate x <= min(y, z) to preserve convexity:
            x <= t && t <= y && t <= z where t is a subunit variable.

            The remainder of the code implements these ideas.
            =#

            # build up the subunit variables
            kstart = k
            for grr in grrs[rid]
                for gid in grr
                    gidx = first(indexin([gid], genes(model)))
                    @constraint(opt_model, t[k] <= y[gidx])
                end
                k += 1
            end
            kend = k - 1

            # total enzyme concentration is sum of minimums
            isozymes = @expression(opt_model, kcats[rid] * sum(t[kstart:kend]))

            if lbs[ridx] >= 0 && ubs[ridx] > 0 # forward only
                @constraint(opt_model, x[ridx] <= isozymes)
                @constraint(opt_model, 0 <= x[ridx])
            elseif lbs[ridx] < 0 && ubs[ridx] <= 0 # reverse only
                @constraint(opt_model, -isozymes <= x[ridx])
                @constraint(opt_model, x[ridx] <= 0)
            elseif lbs[ridx] == ubs[ridx] == 0 # set to zero
                @constraint(opt_model, x[ridx] == 0)
            else # reversible
                @constraint(opt_model, x[ridx] <= isozymes)
                @constraint(opt_model, -isozymes <= x[ridx])
            end
        end
    end
