"""
$(TYPEDSIGNATURES)

A modification that zeroes the bounds of all reactions that would be knocked
out by the combination of specified genes (effectively disabling the
reactions).

A slightly counter-intuitive behavior may occur if knocking out multiple genes:
Because this only changes the reaction bounds, multiple gene knockouts _must_
be specified in a single call to [`knockout`](@ref), because the modifications
have no way to remember which genes are already knocked out and which not.

In turn, having a reaction that can be catalyzed either by Gene1 or by Gene2,
specifying `modifications = [knockout(["Gene1", "Gene2"])]` does indeed disable
the reaction, but `modifications = [knockout("Gene1"), knockout("Gene2")]` does
_not_ disable the reaction (although reactions that depend either only on Gene1
or only on Gene2 are disabled).
"""
knockout(gene_ids::Vector{String}) =
    (model, optmodel) -> _do_knockout(model, optmodel, gene_ids)

"""
$(TYPEDSIGNATURES)

A helper variant of [`knockout`](@ref) for a single gene.
"""
knockout(gene_id::String) = knockout([gene_id])

"""
$(TYPEDSIGNATURES)

Internal helper for knockouts on generic AbstractMetabolicModels. This can be
overloaded so that the knockouts may work differently (more efficiently) with
other models.
"""
function _do_knockout(model::AbstractMetabolicModel, opt_model, gene_ids::Vector{String})
    #TODO this should preferably work on reactions. Make it a wrapper.
    KOs = Set(gene_ids)
    for (ridx, rid) in enumerate(variable_ids(model))
        if eval_reaction_gene_association(model, rid, falses = KOs) == false # also tests for nothing!
            set_optmodel_bound!(ridx, opt_model, lower_bound = 0, upper_bound = 0)
        end
    end
end
