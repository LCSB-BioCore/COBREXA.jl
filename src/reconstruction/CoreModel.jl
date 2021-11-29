"""
    add_reactions!(model::CoreModel, rxns::Vector{Reaction})

Add `rxns` to `model` efficiently. The model must already contain the metabolites used by
`rxns` in the model.
"""
function add_reactions!(model::CoreModel, rxns::Vector{Reaction})
    I = Int64[] # rows
    J = Int64[] # cols
    V = Float64[] # values
    cs = zeros(length(rxns))
    lbs = zeros(length(rxns))
    ubs = zeros(length(rxns))
    for (j, rxn) in enumerate(rxns)
        req = rxn.metabolites
        for (i, v) in zip(indexin(keys(req), metabolites(model)), values(req))
            push!(J, j)
            push!(I, i)
            push!(V, v)
        end
        push!(model.rxns, rxn.id)
        lbs[j] = rxn.lb
        ubs[j] = rxn.ub
        cs[j] = rxn.objective_coefficient
    end
    Sadd = sparse(I, J, V, n_metabolites(model), length(rxns))
    model.S = [model.S Sadd]
    model.c = dropzeros([model.c; cs])
    model.xu = ubs
    model.xl = lbs
    return nothing
end

"""
    add_reaction!(model::CoreModel, rxn::Reaction)

Add `rxn` to `model`. The model must already contain the metabolites used by
`rxn` in the model.
"""
add_reaction!(model::CoreModel, rxn::Reaction) = add_reactions!(model, [rxn])

"""
    add_reactions(
        m::CoreModel,
        s::VecType,
        b::VecType,
        c::AbstractFloat,
        xl::AbstractFloat,
        xu::AbstractFloat;
        check_consistency = false,
    )

Add reaction(s) to a `CoreModel` model `m`.
"""
function add_reactions(
    m::CoreModel,
    s::VecType,
    b::VecType,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat;
    check_consistency = false,
)
    return add_reactions(
        m,
        sparse(reshape(s, (length(s), 1))),
        sparse(b),
        sparse([c]),
        sparse([xl]),
        sparse([xu]),
        check_consistency = check_consistency,
    )
end

"""
    add_reactions(
        m::CoreModel,
        s::VecType,
        b::VecType,
        c::AbstractFloat,
        xl::AbstractFloat,
        xu::AbstractFloat,
        rxn::String,
        mets::K;
        check_consistency = false,
    )
"""
function add_reactions(
    m::CoreModel,
    s::VecType,
    b::VecType,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat,
    rxn::String,
    mets::StringVecType;
    check_consistency = false,
)
    return add_reactions(
        m,
        sparse(reshape(s, (length(s), 1))),
        sparse(b),
        sparse([c]),
        sparse([xl]),
        sparse([xu]),
        [rxn],
        mets,
        check_consistency = check_consistency,
    )
end

"""
    add_reactions(
        m::CoreModel,
        Sp::MatType,
        b::VecType,
        c::VecType,
        xl::VecType,
        xu::VecType;
        check_consistency = false,
    )
"""
function add_reactions(
    m::CoreModel,
    Sp::MatType,
    b::VecType,
    c::VecType,
    xl::VecType,
    xu::VecType;
    check_consistency = false,
)
    rxns = ["r$x" for x = length(m.rxns)+1:length(m.rxns)+length(xu)]
    mets = ["m$x" for x = length(m.mets)+1:length(m.mets)+size(Sp)[1]]
    return add_reactions(
        m,
        Sp,
        b,
        c,
        xl,
        xu,
        rxns,
        mets,
        check_consistency = check_consistency,
    )
end

"""
    add_reactions(m1::CoreModel, m2::CoreModel; check_consistency = false)

Add all reactions from `m2` to `m1`.
"""
function add_reactions(m1::CoreModel, m2::CoreModel; check_consistency = false)
    return add_reactions(
        m1,
        m2.S,
        m2.b,
        m2.c,
        m2.xl,
        m2.xu,
        m2.rxns,
        m2.mets,
        check_consistency = check_consistency,
    )
end

"""
    add_reactions(
        m::CoreModel,
        Sp::MatType,
        b::VecType,
        c::VecType,
        xl::VecType,
        xu::VecType,
        rxns::StringVecType,
        mets::StringVecType;
        check_consistency = false,
    )
"""
function add_reactions(
    m::CoreModel,
    Sp::MatType,
    b::VecType,
    c::VecType,
    xl::VecType,
    xu::VecType,
    rxns::StringVecType,
    mets::StringVecType;
    check_consistency = false,
)
    Sp = sparse(Sp)
    b = sparse(b)
    c = sparse(c)
    xl = collect(xl)
    xu = collect(xu)

    all([length(b), length(mets)] .== size(Sp, 1)) ||
        throw(DimensionMismatch("inconsistent number of metabolites"))
    all(length.([c, xl, xu, rxns]) .== size(Sp, 2)) ||
        throw(DimensionMismatch("inconsistent number of reactions"))

    new_reactions = findall(Bool[!(rxn in m.rxns) for rxn in rxns])
    new_metabolites = findall(Bool[!(met in m.mets) for met in mets])

    if check_consistency
        (newReactions1, newMetabolites1) = verify_consistency(
            m,
            Sp,
            b,
            c,
            xl,
            xu,
            rxns,
            mets,
            new_reactions,
            new_metabolites,
        )
    end

    new_mets = vcat(m.mets, mets[new_metabolites])

    zero_block = spzeros(length(new_metabolites), n_reactions(m))
    ext_s = vcat(sparse(m.S), zero_block)

    mapping = [findfirst(isequal(met), new_mets) for met in mets]
    (I, J, elements) = findnz(sparse(Sp[:, new_reactions]))
    ext_sp = spzeros(length(new_mets), length(new_reactions))
    for (k, i) in enumerate(I)
        new_i = mapping[i]
        ext_sp[new_i, J[k]] = elements[k]
    end

    new_s = hcat(ext_s, ext_sp)
    newb = vcat(m.b, b[new_metabolites])
    newc = vcat(m.c, c[new_reactions])
    newxl = vcat(m.xl, xl[new_reactions])
    newxu = vcat(m.xu, xu[new_reactions])
    new_rxns = vcat(m.rxns, rxns[new_reactions])
    new_lp = CoreModel(new_s, newb, newc, newxl, newxu, new_rxns, new_mets)

    if check_consistency
        return (new_lp, new_reactions, new_metabolites)
    else
        return new_lp
    end
end

"""
    verify_consistency(
        m::CoreModel,
        Sp::M,
        b::V,
        c::V,
        xl::B,
        xu::B,
        names::K,
        mets::K,
        new_reactions,
        new_metabolites,
    ) where {M<:MatType,V<:VecType,B<:VecTypeK<:StringVecType}

Check the consistency of given reactions with existing reactions in `m`.

TODO: work in progress, doesn't return consistency status.
"""
function verify_consistency(
    m::CoreModel,
    Sp::M,
    b::V,
    c::V,
    xl::B,
    xu::B,
    names::K,
    mets::K,
    new_reactions,
    new_metabolites,
) where {M<:MatType,V<:VecType,B<:VecType,K<:StringVecType}

    if !isempty(new_reactions)
        statuses = Vector{ReactionStatus}(undef, length(names))
        for (i, name) in enumerate(names)
            rxn_index = findfirst(isequal(name), m.rxns)
            reaction = Sp[:, i]
            stoich_index = findfirst(Bool[reaction == m.S[:, j] for j = 1:size(m.S, 2)])
            if isnothing(rxn_index) & isnothing(stoich_index)
                statuses[i] = ReactionStatus(false, 0, "new")
            end

            if !isnothing(rxn_index) & isnothing(stoich_index)
                statuses[i] = ReactionStatus(true, 0, "same name")
            end

            if isnothing(rxn_index) & !isnothing(stoich_index)
                statuses[i] = ReactionStatus(true, 0, "same stoichiometry")
            end

            if !isnothing(rxn_index) & !isnothing(stoich_index)
                statuses[i] = ReactionStatus(true, 0, "same name, same stoichiometry")
            end
        end
    end

    return (new_reactions, new_metabolites)
end

@_change_bounds_fn CoreModel Int inplace begin
    isnothing(lower) || (model.xl[rxn_idx] = lower)
    isnothing(upper) || (model.xu[rxn_idx] = upper)
    nothing
end

@_change_bounds_fn CoreModel Int inplace plural begin
    for (i, l, u) in zip(rxn_idxs, lower, upper)
        change_bound!(model, i, lower = l, upper = u)
    end
end

@_change_bounds_fn CoreModel Int begin
    change_bounds(model, [rxn_idx], lower = [lower], upper = [upper])
end

@_change_bounds_fn CoreModel Int plural begin
    n = copy(model)
    n.xl = copy(n.xl)
    n.xu = copy(n.xu)
    change_bounds!(n, rxn_idxs, lower = lower, upper = upper)
    n
end

@_change_bounds_fn CoreModel String inplace begin
    change_bounds!(model, [rxn_id], lower = [lower], upper = [upper])
end

@_change_bounds_fn CoreModel String inplace plural begin
    change_bounds!(
        model,
        Vector{Int}(indexin(rxn_ids, reactions(model))),
        lower = lower,
        upper = upper,
    )
end

@_change_bounds_fn CoreModel String begin
    change_bounds(model, [rxn_id], lower = [lower], upper = [upper])
end

@_change_bounds_fn CoreModel String plural begin
    change_bounds(
        model,
        Int.(indexin(rxn_ids, reactions(model))),
        lower = lower,
        upper = upper,
    )
end

@_remove_fn reaction CoreModel Int inplace begin
    remove_reactions!(model, [reaction_idx])
end

@_remove_fn reaction CoreModel Int inplace plural begin
    mask = .!in.(1:n_reactions(model), Ref(reaction_idxs))
    model.S = model.S[:, mask]
    model.c = model.c[mask]
    model.xl = model.xl[mask]
    model.xu = model.xu[mask]
    model.rxns = model.rxns[mask]
    nothing
end

@_remove_fn reaction CoreModel Int begin
    remove_reactions(model, [reaction_idx])
end

@_remove_fn reaction CoreModel Int plural begin
    n = copy(model)
    remove_reactions!(n, reaction_idxs)
    return n
end

@_remove_fn reaction CoreModel String inplace begin
    remove_reactions!(model, [reaction_id])
end

@_remove_fn reaction CoreModel String inplace plural begin
    remove_reactions!(model, Int.(indexin(reaction_ids, reactions(model))))
end

@_remove_fn reaction CoreModel String begin
    remove_reactions(model, [reaction_id])
end

@_remove_fn reaction CoreModel String plural begin
    remove_reactions(model, Int.(indexin(reaction_ids, reactions(model))))
end

@_remove_fn metabolite CoreModel Int inplace begin
    remove_metabolites!(model, [metabolite_idx])
end

@_remove_fn metabolite CoreModel Int plural inplace begin
    remove_reactions!(
        model,
        [
            ridx for ridx in 1:n_reactions(model) if
            any(in.(findnz(model.S[:, ridx])[1], Ref(metabolite_idxs)))
        ],
    )
    mask = .!in.(1:n_metabolites(model), Ref(metabolite_idxs))
    model.S = model.S[mask, :]
    model.b = model.b[mask]
    model.mets = model.mets[mask]
    nothing
end

@_remove_fn metabolite CoreModel Int begin
    remove_metabolites(model, [metabolite_idx])
end

@_remove_fn metabolite CoreModel Int plural begin
    n = deepcopy(model) #everything gets changed anyway
    remove_metabolites!(n, metabolite_idxs)
    return n
end

@_remove_fn metabolite CoreModel String inplace begin
    remove_metabolites!(model, [metabolite_id])
end

@_remove_fn metabolite CoreModel String inplace plural begin
    remove_metabolites!(model, Int.(indexin(metabolite_ids, metabolites(model))))
end

@_remove_fn metabolite CoreModel String begin
    remove_metabolites(model, [metabolite_id])
end

@_remove_fn metabolite CoreModel String plural begin
    remove_metabolites(model, Int.(indexin(metabolite_ids, metabolites(model))))
end

"""
    change_objective!(
        model::CoreModel,
        rxn_idxs::Vector{Int};
        weights = ones(length(rxn_idxs)),
    )

Change the objective for `model` to reaction(s) with indices `rxn_idxs`, optionally
specifying their `weights`. By default, assume equal weights. If no objective exists in
model, sets objective. Note, also accepts `String` or `Vector{String}` of reaction ids.
"""
function change_objective!(
    model::CoreModel,
    rxn_idxs::Vector{Int};
    weights = ones(length(rxn_idxs)),
)
    model.c = spzeros(length(model.c))
    model.c[rxn_idxs] .= weights
    nothing
end

change_objective!(model::CoreModel, rxn_idx::Int) = change_objective!(model, [rxn_idx])

function change_objective!(
    model::CoreModel,
    rxn_ids::Vector{String};
    weights = ones(length(rxn_ids)),
)
    idxs = indexin(rxn_ids, reactions(model))
    any(isnothing(idx) for idx in idxs) &&
        throw(DomainError(rxn_ids, "Some reaction ids were not found in model."))
    change_objective!(model, Int.(idxs); weights)
end

change_objective!(model::CoreModel, rxn_id::String) = change_objective!(model, [rxn_id])
