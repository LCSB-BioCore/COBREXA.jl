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
    remove_reactions(m::CoreModel, rxns::Vector{Int})

Remove reaction(s) from a `CoreModel`.

Also removes any metabolites not involved in any reaction after the deletion.

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

"""
    change_bounds!(
        model::CoreModel,
        rxns::Vector{Int};
        xl::V = Float64[],
        xu::V = Float64[],
    ) where {V<:VecType}

Change the lower and/or upper bounds ('xl' and 'xu') for given reactions.

"""
function change_bounds!(
    model::CoreModel,
    rxns::Vector{Int};
    xl::V = Float64[],
    xu::V = Float64[],
) where {V<:VecType}
    found = [index ∈ 1:n_reactions(model) for index in rxns]
    length(rxns[found]) == length(unique(rxns[found])) ||
        error("`rxns` appears to contain duplicates")

    if !isempty(xl)
        length(rxns) == length(xl) ||
            throw(DimensionMismatch("`rxns` size doesn't match with `xl`"))
        model.xl[rxns[found]] = xl[found]
    end
    if !isempty(xu)
        length(rxns) == length(xu) ||
            throw(DimensionMismatch("`rxns` size doesn't match with `xu`"))
        model.xu[rxns[found]] = xu[found]
    end
end

"""
    change_bounds!(
        model::CoreModel,
        rxns::Vector{String};
        xl::V = Float64[],
        xu::V = Float64[],
    ) where {V<:VecType}

"""
function change_bounds!(
    model::CoreModel,
    rxns::Vector{String};
    xl::V = Float64[],
    xu::V = Float64[],
) where {V<:VecType}
    found = [name ∈ model.rxns for name in rxns]
    rxn_indices = zeros(Int, length(rxns))
    rxn_indices[found] = [findfirst(isequal(name), model.rxns) for name in rxns[found]]
    change_bounds!(model, rxn_indices, xl = xl, xu = xu)
end
