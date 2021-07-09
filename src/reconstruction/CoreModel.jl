"""
    add_reactions(
        m::CoreModel,
        s::V1,
        b::V2,
        c::AbstractFloat,
        xl::AbstractFloat,
        xu::AbstractFloat;
        check_consistency = false,
    ) where {V1<:VecType,V2<:VecType}

Add reaction(s) to a `CoreModel` model `m`.

"""
function add_reactions(
    m::CoreModel,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType}
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
        s::V1,
        b::V2,
        c::AbstractFloat,
        xl::AbstractFloat,
        xu::AbstractFloat,
        rxn::String,
        mets::K;
        check_consistency = false,
    ) where {V1<:VecType,V2<:VecType,K<:StringVecType}

"""
function add_reactions(
    m::CoreModel,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat,
    rxn::String,
    mets::K;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType,K<:StringVecType}
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
        Sp::M,
        b::V,
        c::V,
        xl::V,
        xu::V;
        check_consistency = false,
    ) where {M<:MatType,V<:VecType}

"""
function add_reactions(
    m::CoreModel,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V;
    check_consistency = false,
) where {M<:MatType,V<:VecType}
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
        Sp::M,
        b::V,
        c::V,
        xl::V,
        xu::V,
        rxns::K,
        mets::K;
        check_consistency = false,
    ) where {M<:MatType,V<:VecType,K<:StringVecType}

"""
function add_reactions(
    m::CoreModel,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    rxns::K,
    mets::K;
    check_consistency = false,
) where {M<:MatType,V<:VecType,K<:StringVecType}

    Sp = sparse(Sp)
    b = sparse(b)
    c = sparse(c)
    xl = sparse(xl)
    xu = sparse(xu)

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
    #new_c = hcat(m.C, spzeros(size(m.C, 1), length(new_reactions)))
    newc = vcat(m.c, c[new_reactions])
    newxl = vcat(m.xl, xl[new_reactions])
    newxu = vcat(m.xu, xu[new_reactions])
    new_rxns = vcat(m.rxns, rxns[new_reactions])
    #new_lp = CoreModel(new_s, newb, new_c, m.cl, m.cu, newc, newxl, newxu, new_rxns, new_mets)
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
        xl::V,
        xu::V,
        names::K,
        mets::K,
        new_reactions,
        new_metabolites,
    ) where {M<:MatType,V<:VecType,K<:StringVecType}

Check the consistency of given reactions with existing reactions in `m`.

TODO: work in progress, doesn't return consistency status.

"""
function verify_consistency(
    m::CoreModel,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    names::K,
    mets::K,
    new_reactions,
    new_metabolites,
) where {M<:MatType,V<:VecType,K<:StringVecType}

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

"""
    remove_metabolites(model::CoreModel, metabolites)

Removes a set of `metabolites` from the `model` of type `CoreModel` and returns
a new `CoreModel` without those metabolites. Here, `metabolites` can be either a
string, a vector of strings, an index or a vector of indices. Also removes any
reactions that have no associated metabolites after the metabolites have been
removed.

# Example
```
model = load_model(CoreModel, "e_coli_core.json")

m1 = remove_metabolites(model, ["glc__D_e", "for_c"])
m2 = remove_metabolites(model, "glc__D_e")
m3 = remove_metabolites(model, indexin(["glc__D_e", "for_c"], metabolites(model)))
m4 = remove_metabolites(model, first(indexin(["glc__D_e"], metabolites(model))))
```
"""
function remove_metabolites(model::CoreModel, mets)
    mets_to_keep = filter(x -> x ∉ mets, 1:n_metabolites(model))
    temp_S = model.S[mets_to_keep, :]

    (I, rxns_to_keep, val) = findnz(temp_S)
    sort!(rxns_to_keep)
    unique!(rxns_to_keep)
    new_S = model.S[mets_to_keep, rxns_to_keep]
    new_b = model.b[mets_to_keep]
    new_c = model.c[rxns_to_keep]
    new_lbs = model.xl[rxns_to_keep]
    new_ubs = model.xu[rxns_to_keep]
    new_rxns = model.rxns[rxns_to_keep]
    new_mets = model.mets[mets_to_keep]

    return CoreModel(new_S, new_b, new_c, new_lbs, new_ubs, new_rxns, new_mets)
end

function remove_metabolites(model::CoreModel, met::Integer)
    return remove_metabolites(model, [met])
end

function remove_metabolites(model::CoreModel, met::String)
    return remove_metabolites(model, [met])
end

function remove_metabolites(model::CoreModel, mets::Vector{String})
    met_indices = filter(!isnothing, indexin(mets, metabolites(model)))
    if isempty(met_indices)
        return model
    else
        return remove_metabolites(model, met_indices)
    end
end

"""
    remove_reactions(m::CoreModel, rxns::Vector{Int})

Removes a set of reactions from a CoreModel.
Also removes the metabolites not involved in any reaction.
"""
function remove_reactions(m::CoreModel, rxns::Vector{Int})
    rxns_to_keep = filter(e -> e ∉ rxns, 1:n_reactions(m))
    temp_s = m.S[:, rxns_to_keep]

    (mets_to_keep, J, val) = findnz(temp_s)
    sort!(mets_to_keep)
    unique!(mets_to_keep)
    new_s = m.S[mets_to_keep, rxns_to_keep]
    newb = m.b[mets_to_keep]
    newc = m.c[rxns_to_keep]
    newxl = m.xl[rxns_to_keep]
    newxu = m.xu[rxns_to_keep]
    new_rxns = m.rxns[rxns_to_keep]
    new_mets = m.mets[mets_to_keep]
    new_model = CoreModel(new_s, newb, newc, newxl, newxu, new_rxns, new_mets)
    return new_model
end

"""
    remove_reactions(m::CoreModel, rxn::Integer)

"""
function remove_reactions(m::CoreModel, rxn::Integer)
    return remove_reactions(m, [rxn])
end

"""
    remove_reactions(m::CoreModel, rxn::String)

"""
function remove_reactions(m::CoreModel, rxn::String)
    return remove_reactions(m, [rxn])
end

"""
    remove_reactions(m::CoreModel, rxns::Vector{String})

"""
function remove_reactions(m::CoreModel, rxns::Vector{String})
    rxn_indices = [findfirst(isequal(name), m.rxns) for name in intersect(rxns, m.rxns)]
    if isempty(rxn_indices)
        return m
    else
        return remove_reactions(m, rxn_indices)
    end
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
        Vector{Int}(indexin(rxn_ids, reactions(model))),
        lower = lower,
        upper = upper,
    )
end
