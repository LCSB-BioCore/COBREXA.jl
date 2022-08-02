"""
    add_knockout_constraints(gene_ids::Vector{String})

A modification that zeroes the bounds of all reactions that would be knocked
out by the combination of specified genes (effectively disabling the
reactions).

A slightly counter-intuitive behavior may occur if knocking out multiple genes:
Because this only changes the reaction bounds, multiple gene knockouts _must_
be specified in a single call to [`add_knockout_constraints`](@ref), because the modifications
have no way to remember which genes are already knocked out and which not.

In turn, having a reaction that can be catalyzed either by Gene1 or by Gene2,
specifying `modifications = [add_knockout_constraints(["Gene1", "Gene2"])]` does indeed disable
the reaction, but `modifications = [add_knockout_constraints("Gene1"), add_knockout_constraints("Gene2")]` does
_not_ disable the reaction (although reactions that depend either only on Gene1
or only on Gene2 are disabled).
"""
add_knockout_constraints(gene_ids::Vector{String}) =
    (model, optmodel) -> _do_knockout_constraints(model, optmodel, gene_ids)

"""
    add_knockout_constraints(gene_id::String)

A helper variant of [`add_knockout_constraints`](@ref) for a single gene.
"""
add_knockout_constraints(gene_id::String) = add_knockout_constraints([gene_id])

"""
    _do_knockout(model::MetabolicModel, opt_model)

Internal helper for knockouts on generic MetabolicModels. This can be
overloaded so that the knockouts may work differently (perhaps more
efficiently) with other model types.
"""
function _do_knockout_constraints(
    model::MetabolicModel,
    opt_model,
    gene_ids::Vector{String},
)
    for (rxn_num, rxn_id) in enumerate(reactions(model))
        rga = reaction_gene_association(model, rxn_id)
        if !isnothing(rga) &&
           all([any(in.(gene_ids, Ref(conjunction))) for conjunction in rga])
            set_optmodel_bound!(rxn_num, opt_model, ub = 0, lb = 0)
        end
    end
end
