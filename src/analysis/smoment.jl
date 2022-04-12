"""
    make_smomentmodel(
        model::StandardModel;
        rid_isozymes = Dict{String, Vector{Isozyme}}(),
    )

Construct an `SMomentModel` model using `model` and `rid_isozymes`.
"""
function make_smomentmodel(
    model::StandardModel;
    rid_isozymes = Dict{String,Vector{Isozyme}}(),
    enzyme_capacity = 0.0,
)

    # check that input data is in correct format for smoment
    if any(length(v) > 1 for v in values(rid_isozymes))
        @warn(
            "For SMOMENT to work correctly, no isozymes are allowed. Call `remove_slow_isozymes!` to fix the input data."
        )
    end

    irrevS, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model)

    #: size of resultant model
    num_reactions = size(irrevS, 2)
    num_metabolites = size(irrevS, 1)
    num_vars = num_reactions + 1

    #: equality lhs
    Se = zeros(1, num_reactions)

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(rid_isozymes, original_rid) && continue
        # these entries have kcats, only one GRR by assumption
        isozyme = first(rid_isozymes[original_rid])
        mw = sum([model.genes[gid].molar_mass * ps for (gid, ps) in isozyme.stoichiometry])
        kcat = contains(rid, "§FOR") ? isozyme.kcat_forward : isozyme.kcat_reverse
        Se[1, col_idx] = -mw / kcat
    end

    S = [
        irrevS zeros(num_metabolites, 1)
        Se 1.0
    ]

    #: equality rhs
    b = zeros(num_metabolites + 1)

    #: find objective
    obj_idx_orig = first(findnz(objective(model))[1])
    obj_id_orig = reactions(model)[obj_idx_orig]
    obj_id = obj_id_orig * "§FOR" # assume forward reaction is objective
    c = spzeros(num_vars)
    obj_idx = reaction_map[obj_id]
    c[obj_idx] = 1.0

    #: bounds
    xl = sparse([lb_fluxes; 0.0])
    xu = sparse([ub_fluxes; enzyme_capacity])

    return SMomentModel(
        reactions(model),
        _order_id_to_idx_dict(reaction_map),
        _order_id_to_idx_dict(metabolite_map),
        c,
        S,
        b,
        xl,
        xu,
    )
end
