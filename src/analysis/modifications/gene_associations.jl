
"""
    change_eflux_constraints(
        args...;
        availability_func = enzyme_availability_eflux,
        reaction_gene_association_func = reaction_gene_association,
        normalize = true,
        kwargs...,
    )

Analysis modification that changes the existing reaction bounds accordingly to
the E-Flux method. Arguments are forwarded to the `availability_func`, by
default [`enzyme_availability_eflux`](@ref), which is used to find the relative
adjustments of the reactions. Reactions with no associated genes are not
modified (the data is retrieved using the `reaction_gene_association_func`, by
default the [`reaction_gene_association`](@ref)).

If `normalize` is true, all relative adjustments are normalized to maximum
1.

For details see: Colijn, Caroline, et al. "Interpreting expression data with
metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid
production." PLoS computational biology 5.8 (2009): e1000489.
"""
change_eflux_constraints(
    args...;
    availability_func = enzyme_availability_eflux,
    reaction_gene_association_func = reaction_gene_association,
    normalize = true,
    kwargs...,
) =
    (model, opt_model) -> begin

        adjustments = Pair{Int,Float64}[]
        for (ridx, rid) in enumerate(reactions(model))
            grr = reaction_gene_association(model, rid)
            isnothing(grr) && continue
            push!(adjustments, rid => availability_function(grr, args...; kwargs...))
        end

        if normalize && !isempty(adjustments)
            max = maximum(adj for (_, a) in adjustments)
            for (_, adj) in adjustments
                adj /= max
            end
        end

        lbs, ubs = get_optmodel_bounds(opt_model)
        for (ridx, adj) in adjustments
            set_optmodel_bound(a, opt_model, lbs[ridx] * adj, ubs[ridx] * adj)
        end
    end
