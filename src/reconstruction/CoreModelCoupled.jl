"""
Adds reactions to the coupled model `m`
"""
function add_reactions(
    m::CoreModelCoupled,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType}
    new_lm = add_reactions(
                m.lm,
                s,
                b,
                c,
                xl,
                xu,
                check_consistency = check_consistency,
            )
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm)-n_reactions(m.lm))),
        m.cl,
        m.cu
    )
end

function add_reactions(
    m::CoreModelCoupled,
    s::V1,
    b::V2,
    c::AbstractFloat,
    xl::AbstractFloat,
    xu::AbstractFloat,
    rxn::String,
    mets::K;
    check_consistency = false,
) where {V1<:VecType,V2<:VecType,K<:StringVecType}
    new_lm = add_reactions(
                m.lm,
                s,
                b,
                c,
                xl,
                xu,
                rxn,
                mets,
                check_consistency = check_consistency
            )
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm)-n_reactions(m.lm))),
        m.cl,
        m.cu
    )
end

function add_reactions(
    m::CoreModelCoupled,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V;
    check_consistency = false,
) where {M<:MatType,V<:VecType}
    new_lm = add_reactions(
                m.lm,
                Sp,
                b,
                c,
                xl,
                xu,
                check_consistency = check_consistency,
            )
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm)-n_reactions(m.lm))),
        m.cl,
        m.cu
    )
end

function add_reactions(m1::CoreModelCoupled, m2::CoreModel; check_consistency = false)
    new_lm = add_reactions(
                m1.lm,
                m2,
                check_consistency=check_consistency
            )
    return CoreModelCoupled(
        new_lm,
        hcat(m1.C, spzeros(size(m1.C, 1), n_reactions(new_lm)-n_reactions(m1.lm))),
        m1.cl,
        m1.cu
    )
end

function add_reactions(
    m::CoreModelCoupled,
    Sp::M,
    b::V,
    c::V,
    xl::V,
    xu::V,
    rxns::K,
    mets::K;
    check_consistency = false,
) where {M<:MatType,V<:VecType,K<:StringVecType}
    new_lm = add_reactions(
            m.lm,
            Sp,
            b,
            c,
            xl,
            xu,
            rxns,
            mets,
            check_consistency = check_consistency
        )
    return CoreModelCoupled(
        new_lm,
        hcat(m.C, spzeros(size(m.C, 1), n_reactions(new_lm)-n_reactions(m.lm))),
        m.cl,
        m.cu
    )
end

"""
Remove a set of reactions from a `CoreModelCoupled`.

Also removes the metabolites not involved in any reaction.
"""
function remove_reactions(m::CoreModelCoupled, rxns::Vector{Int})
    return CoreModelCoupled(
                remove_reactions(m.lm, rxns),
                m.C[:, filter(e -> e ∉ rxns, 1:n_reactions(m))],
                m.cl,
                m.cu
            )
end

function remove_reactions(m::CoreModelCoupled, rxn::Integer)
    return remove_reactions(m, [rxn])
end

function remove_reactions(m::CoreModelCoupled, rxn::String)
    return remove_reactions(m, [rxn])
end

function remove_reactions(m::CoreModelCoupled, rxns::Vector{String})
    rxn_indices = [findfirst(isequal(name), m.lm.rxns) for name in intersect(rxns, m.lm.rxns)]
    if isempty(rxn_indices)
        return m
    else
        return remove_reactions(m, rxn_indices)
    end
end

"""
Add constraints of the following form to a CoreModelCoupled and return a modified one.

The arguments are same as for in-place `add_coupling_constraints!`.
"""
function add_coupling_constraints(m::CoreModelCoupled, args...)
    new_lp = deepcopy(m)
    add_coupling_constraints!(new_lp, args...)
    return new_lp
end

"""
Add constraints to a plain `CoreModel` (converts it to the coupled type)
"""
add_coupling_constraints(m::CoreModel, args...) =
    add_coupling_constraints(convert(CoreModelCoupled, m), args...)

"""
In-place add coupling constraints in form
```
    cₗ ≤ C x ≤ cᵤ
```
"""
function add_coupling_constraints!(
    m::CoreModelCoupled,
    c::V,
    cl::AbstractFloat,
    cu::AbstractFloat,
) where {V<:VecType}
    return add_coupling_constraints!(
        m,
        sparse(reshape(c, (1, length(c)))),
        sparse([cl]),
        sparse([cu]),
    )
end


function add_coupling_constraints!(
    m::CoreModelCoupled,
    C::M,
    cl::V,
    cu::V,
) where {M<:MatType,V<:VecType}

    all([length(cu), length(cl)] .== size(C, 1)) ||
        throw(DimensionMismatch("mismatched numbers of constraints"))
    size(C, 2) == n_reactions(m) ||
        throw(DimensionMismatch("mismatched number of reactions"))

    m.C = vcat(m.C, sparse(C))
    m.cl = vcat(m.cl, sparse(cl))
    m.cu = vcat(m.cu, sparse(cu))
end


"""
Remove coupling constraints from the linear model and return the modified model.

Arguments are the same as for in-place version `remove_coupling_constraints!`.
"""
function remove_coupling_constraints(m::CoreModelCoupled, args...)
    new_model = deepcopy(m)
    remove_coupling_constraints!(new_model, args...)
    return new_model
end


"""
Removes a set of coupling constraints from a CoreModelCoupled in-place.
"""
function remove_coupling_constraints!(m::CoreModelCoupled, constraint::Int)
    remove_coupling_constraints!(m, [constraint])
end


function remove_coupling_constraints!(m::CoreModelCoupled, constraints::Vector{Int})
    to_be_kept = filter(e -> e ∉ constraints, 1:n_coupling_constraints(m))
    m.C = m.C[to_be_kept, :]
    m.cl = m.cl[to_be_kept]
    m.cu = m.cu[to_be_kept]
end


"""
Change the lower and/or upper bounds ('cl' and 'cu') for given coupling constraints
"""
function change_coupling_bounds!(
    model::CoreModelCoupled,
    constraints::Vector{Int};
    cl::V = Float64[],
    cu::V = Float64[],
) where {V<:VecType}
    found = [index ∈ 1:n_coupling_constraints(model) for index in constraints]
    red_constraints = constraints[found]

    length(red_constraints) == length(unique(red_constraints)) ||
        error("`constraints` appears to contain duplicates")
    if !isempty(cl)
        length(constraints) == length(cl) ||
            throw(DimensionMismatch("`constraints` size doesn't match with `cl`"))
        model.cl[red_constraints] = cl[found]
    end

    if !isempty(cu)
        length(constraints) == length(cu) ||
            throw(DimensionMismatch("`constraints` size doesn't match with `cu`"))
        model.cu[red_constraints] = cu[found]
    end
end

"""
Return indices of exchange reactions.

Exchange reactions are identified based on most commonly used prefixes.
"""
function find_exchange_reactions(
    model::CoreModelCoupled;
    exclude_biomass = false,
    biomass_str::String = "biomass",
    exc_prefs = ["EX_"; "Exch_"; "Ex_"],
)
    return find_exchange_reactions(
        model.lm;
        exclude_biomass = exclude_biomass,
        biomass_str = biomass_str,
        exc_prefs = exc_prefs,
    )
end

"""
Return indices of exchange metabolites.

Exchange metabolites are identified based on the exchange reactions return
by `find_exchange_reactions`.
"""
function find_exchange_metabolites(
    model::CoreModelCoupled;
    exclude_biomass = false,
    biomass_str::String = "biomass",
    exc_prefs = ["EX_"; "Exch_"; "Ex_"],
)
    return find_exchange_metabolites(
        model.lm;
        exclude_biomass = exclude_biomass,
        biomass_str = biomass_str,
        exc_prefs = exc_prefs,
    )
end

"""
Change the lower and/or upper bounds ('xl' and 'xu') for given reactions
"""
function change_bounds!(
    model::CoreModelCoupled,
    rxns::Vector{Int};
    xl::V = Float64[],
    xu::V = Float64[],
) where {V<:VecType}
    change_bounds!(
        model.lm,
        rxns,
        xl = xl,
        xu = xu
    )
end

function change_bounds!(
    model::CoreModelCoupled,
    rxns::Vector{String};
    xl::V = Float64[],
    xu::V = Float64[],
) where {V<:VecType}
    change_bounds!(
        model.lm,
        rxns,
        xl = xl,
        xu = xu
    )
end
