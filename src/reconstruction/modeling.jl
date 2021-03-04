"""
Adds reactions to the model `m`
"""
function addReactions(
    m::LinearModel,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat;
    checkConsistency = false,
) where {V1<:VT,V2<:VT}
    return addReactions(
        m,
        sparse(reshape(s, (length(s), 1))),
        sparse(b),
        sparse([c]),
        sparse([xl]),
        sparse([xu]),
        checkConsistency = checkConsistency,
    )
end


function addReactions(
    m::LinearModel,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat,
    rxn::String,
    mets::K;
    checkConsistency = false,
) where {V1<:VT,V2<:VT,K<:ST}
    return addReactions(
        m,
        sparse(reshape(s, (length(s), 1))),
        sparse(b),
        sparse([c]),
        sparse([xl]),
        sparse([xu]),
        [rxn],
        mets,
        checkConsistency = checkConsistency,
    )
end

function addReactions(
    m::LinearModel,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V;
    checkConsistency = false,
) where {M<:MT,V<:VT}
    rxns = ["r$x" for x = length(m.rxns)+1:length(m.rxns)+length(xu)]
    mets = ["m$x" for x = length(m.mets)+1:length(m.mets)+size(Sp)[1]]
    return addReactions(
        m,
        Sp,
        b,
        c,
        xl,
        xu,
        rxns,
        mets,
        checkConsistency = checkConsistency,
    )
end


function addReactions(m1::LinearModel, m2::LinearModel; checkConsistency = false)
    return addReactions(
        m1,
        m2.S,
        m2.b,
        m2.c,
        m2.xl,
        m2.xu,
        m2.rxns,
        m2.mets,
        checkConsistency = checkConsistency,
    )
end


function addReactions(
    m::LinearModel,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    rxns::K,
    mets::K;
    checkConsistency = false,
) where {M<:MT,V<:VT,K<:ST}

    Sp = sparse(Sp)
    b = sparse(b)
    c = sparse(c)
    xl = sparse(xl)
    xu = sparse(xu)

    checkInputDimensions(Sp, b, c, xl, xu, rxns, mets)

    newReactions = findall(Bool[!(rxn in m.rxns) for rxn in rxns])
    newMetabolites = findall(Bool[!(met in m.mets) for met in mets])


    if checkConsistency
        (newReactions1, newMetabolites1) =
            verifyConsistency(m, Sp, b, c, xl, xu, rxns, mets, newReactions, newMetabolites)
    end

    newMets = vcat(m.mets, mets[newMetabolites])

    zeroBlock = spzeros(length(newMetabolites), nReactions(m))
    extS = vcat(sparse(m.S), zeroBlock)

    mapping = [findfirst(isequal(met), newMets) for met in mets]
    (I, J, elements) = findnz(sparse(Sp[:, newReactions]))
    extSp = spzeros(length(newMets), length(newReactions))
    for (k, i) in enumerate(I)
        newI = mapping[i]
        extSp[newI, J[k]] = elements[k]
    end

    newS = hcat(extS, extSp)
    newb = vcat(m.b, b[newMetabolites])
    newC = hcat(m.C, spzeros(size(m.C, 1), length(newReactions)))
    newc = vcat(m.c, c[newReactions])
    newxl = vcat(m.xl, xl[newReactions])
    newxu = vcat(m.xu, xu[newReactions])
    newRxns = vcat(m.rxns, rxns[newReactions])
    newLp = LinearModel(newS, newb, newC, m.cl, m.cu, newc, newxl, newxu, newRxns, newMets)

    if checkConsistency
        return (newLp, newReactions, newMetabolites)
    else
        return newLp
    end
end

"""
Verifies that vectors and matrices have the expected dimensions.
"""
function checkInputDimensions(
    S::M1,
    b::V,
    C::M2,
    cl::V,
    cu::V,
    c::V,
    xl::V,
    xu::V,
    rxns::K,
    mets::K,
) where {M1<:MT,M2<:MT,V<:VT,K<:ST}
    n_c = length(c)

    length(cu) == length(cl) ||
        throw(DimensionMismatch("`cl` and `cu` don't have the same size"))
    size(C) == (length(cu), n_c) ||
        throw(DimensionMismatch("`C` shape doesn't match with `cu` and `c`"))

    checkInputDimensions(S, b, c, xl, xu, rxns, mets)

end

function checkInputDimensions(
    S::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    rxns::K,
    mets::K,
) where {M<:MT,V<:VT,K<:ST}

    n_c = length(c)
    n_b = length(b)

    length(xu) == length(xl) ||
        throw(DimensionMismatch("`xl` and `xu` don't have the same size"))
    n_c == length(xl) || throw(DimensionMismatch("`c` doesn't have the same size as `xl`"))

    size(S) == (n_b, n_c) ||
        throw(DimensionMismatch("`S` shape doesn't match with `c` and `mets`"))

    length(rxns) == n_c || throw(DimensionMismatch("`rxns` size doesn't match with `S`"))
    length(mets) == n_b || throw(DimensionMismatch("`mets` size doesn't match with `S`"))
end

"""
Verifies the consistency of a given model
"""
function verifyConsistency(
    m::LinearModel,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    names::K,
    mets::K,
    newReactions,
    newMetabolites,
) where {M<:MT,V<:VT,K<:ST}

    if !isempty(newReactions)
        statuses = Array{ReactionStatus}(undef, length(names))
        for (i, name) in enumerate(names)
            rxnIndex = findfirst(isequal(name), m.rxns)
            reaction = Sp[:, i]
            stoichIndex = findfirst(Bool[reaction == m.S[:, j] for j = 1:size(m.S, 2)])
            if isnothing(rxnIndex) & isnothing(stoichIndex)
                statuses[i] = ReactionStatus(false, 0, "new")
            end

            if !isnothing(rxnIndex) & isnothing(stoichIndex)
                statuses[i] = ReactionStatus(true, 0, "same name")
            end

            if isnothing(rxnIndex) & !isnothing(stoichIndex)
                statuses[i] = ReactionStatus(true, 0, "same stoichiometry")
            end

            if !isnothing(rxnIndex) & !isnothing(stoichIndex)
                statuses[i] = ReactionStatus(true, 0, "same name, same stoichiometry")
            end
        end
    end

    return (newReactions, newMetabolites)
end

"""
Removes a set of reactions from a LinearModel.
Also removes the metabolites not involved in any reaction.
"""
function removeReactions(m::LinearModel, rxns::Vector{Int})
    rxnsToKeep = filter(e -> e ∉ rxns, 1:nReactions(m))
    tempS = m.S[:, rxnsToKeep]

    (metsToKeep, J, val) = findnz(tempS)
    sort!(metsToKeep)
    unique!(metsToKeep)
    newS = m.S[metsToKeep, rxnsToKeep]
    newb = m.b[metsToKeep]
    newC = m.C[:, rxnsToKeep]
    newc = m.c[rxnsToKeep]
    newxl = m.xl[rxnsToKeep]
    newxu = m.xu[rxnsToKeep]
    newRxns = m.rxns[rxnsToKeep]
    newMets = m.mets[metsToKeep]
    newModel =
        LinearModel(newS, newb, newC, m.cl, m.cu, newc, newxl, newxu, newRxns, newMets)
    return newModel
end

function removeReactions(m::LinearModel, rxn::Integer)
    return removeReactions(m, [rxn])
end


function removeReactions(m::LinearModel, rxn::String)
    return removeReactions(m, [rxn])
end


function removeReactions(m::LinearModel, rxns::Array{String,1})
    rxnIndices = [findfirst(isequal(name), m.rxns) for name in intersect(rxns, m.rxns)]
    if isempty(rxnIndices)
        return m
    else
        return removeReactions(m, rxnIndices)
    end
end


"""
Returns indices of exchange reactions.
Exchange reactions are identified based on most commonly used prefixes.
"""
function findExchangeReactions(
    model::LinearModel;
    excludeBiomass = false,
    biomassStr::String = "biomass",
    excPrefs = ["EX_"; "Exch_"; "Ex_"],
)
    isExc = falses(nReactions(model))
    for pref in excPrefs
        isExc = isExc .| startswith.(model.rxns, pref)
    end
    excInds = findall(isExc)
    if excludeBiomass
        biomInds = findall(x -> occursin(biomassStr, x), model.rxns)
        excInds = setdiff(excInds, biomInds)
    end
    return excInds
end

"""
Returns indices of exchanged metabolites, ie, the outermost metabolites in the network
In practice returns the metabolites consumed by the reactions given by `findExchangeReactions`
and if called with the same arguments, the two outputs correspond.
"""
function findExchangeMetabolites(
    model::LinearModel;
    excludeBiomass = false,
    biomassStr::String = "biomass",
    excPrefs = ["EX_"; "Exch_"; "Ex_"],
)
    excRxnInds = findExchangeReactions(
        model,
        excludeBiomass = excludeBiomass,
        biomassStr = biomassStr,
        excPrefs = excPrefs,
    )
    excMetInds = [findfirst(x -> x == -1, model.S[:, j]) for j in excRxnInds]
    return excMetInds
end


"""
Change the lower and/or upper bounds ('xl' and 'xu') for given reactions
"""
function changeBounds!(
    model::LinearModel,
    rxns::Vector{Int};
    xl::V = Array{Float64}(undef, 0),
    xu::V = Array{Float64}(undef, 0),
) where {V<:VT}
    found = [index ∈ 1:nReactions(model) for index in rxns]
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


function changeBounds!(
    model::LinearModel,
    rxns::Array{String,1};
    xl::V = Array{Float64}(undef, 0),
    xu::V = Array{Float64}(undef, 0),
) where {V<:VT}
    found = [name ∈ model.rxns for name in rxns]
    rxnIndices = zeros(Int, length(rxns))
    rxnIndices[found] = [findfirst(isequal(name), model.rxns) for name in rxns[found]]
    changeBounds!(model, rxnIndices, xl = xl, xu = xu)
end
