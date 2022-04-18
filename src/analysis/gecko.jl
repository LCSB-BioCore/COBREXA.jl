"""
    make_gecko_model(
        model::MetabolicModel;
        reaction_isozymes::Union{Function,Dict{String,Vector{Isozyme}}}
        gene_product_bounds::Union{Function,Dict{String,Tuple{Float64,Float64}}},
        gene_product_molar_mass::Union{Function,Dict{String,Float64}},
        gene_mass_group::Union{Function,Dict{String,String}} = _ -> "uncategorized",
        gene_mass_group_bound::Union{Function,Dict{String,Float64}},
        relaxed_arm_reaction_bounds = false,
    )

Wrap a model into a [`GeckoModel`](@ref), following the structure given by
GECKO algorithm (see [`GeckoModel`](@ref) documentation for details).

# Arguments

- `reaction_isozymes` is a function that returns a vector of [`Isozyme`](@ref)s
  for each reaction, or empty vector if the reaction is not enzymatic.
- `gene_product_bounds` is a function that returns lower and upper bound for
  concentration for a given gene product (specified by the same string gene ID as in
  `reaction_isozymes`), as `Tuple{Float64,Float64}`.
- `gene_product_molar_mass` is a function that returns a numeric molar mass of
  a given gene product specified by string gene ID.
- `gene_mass_group` is a function that returns a string group identifier for a
  given gene product, again specified by string gene ID. By default, all gene
  products belong to group `"uncategorized"` which is the behavior of original
  GECKO.
- `gene_mass_group_bound` is a function that returns the maximum mass for a given
  mass group.
- `relaxed_arm_reaction_bounds` is a boolean flag that relaxes the constraints
  on the "arm" reactions specified by GECKO. By default (value `false`), there
  is a separate constraint that limits the total flux through forward-direction
  reaction for all its isozymes (ensuring that the sum of forward rates is less
  than "global" upper bound), and another separate constraint that limits the
  total flux through reverse-direction reaction isozymes. Value `true` groups
  both forward and reverse reactions in a single constraint, allowing the total
  forward flux to be actually greater than the upper bound IF the reverse flux
  can balance it to fit into the upper and lower bound constraints (in turn,
  more enzyme can be "wasted" by a reaction that runs in both directions).

Alternatively, all function arguments may be also supplied as dictionaries that
provide the same data lookup.
"""
function make_gecko_model(
    model::MetabolicModel;
    reaction_isozymes::Union{Function,Dict{String,Vector{Isozyme}}},
    gene_product_bounds::Union{Function,Dict{String,Tuple{Float64,Float64}}},
    gene_product_molar_mass::Union{Function,Dict{String,Float64}},
    gene_mass_group::Union{Function,Dict{String,String}} = _ -> "uncategorized",
    gene_mass_group_bound::Union{Function,Dict{String,Float64}},
    relaxed_arm_reaction_bounds = false,
)
    ris_ =
        reaction_isozymes isa Function ? reaction_isozymes : (rid -> get(reaction_isozymes, rid, []))
    gpb_ =
        gene_product_bounds isa Function ? gene_product_bounds :
        (gid -> gene_product_bounds[gid])
    gpmm_ =
        gene_product_molar_mass isa Function ? gene_product_molar_mass :
        (gid -> gene_product_molar_mass[gid])
    gmg_ = gene_mass_group isa Function ? gene_mass_group : (gid -> gene_mass_group[gid])
    gmgb_ = gene_mass_group_bound isa Function ? gene_mass_group_bound : (grp -> gene_mass_group_bound[grp])
    # ...it would be nicer to have an overload for this, but kwargs can't be used for dispatch

    columns = Vector{_gecko_column}()
    coupling_row_reaction = Int[]
    coupling_row_gene_product = Int[]
    coupling_row_mass_group = String[]

    gids = genes(model)
    (lbs, ubs) = bounds(model)
    rids = reactions(model)

    gene_name_lookup = Dict(gids .=> 1:length(gids))
    gene_row_lookup = Dict{Int,Int}()
    mass_group_lookup = Dict{String,Int}()

    for i = 1:n_reactions(model)
        isozymes = ris_(rids[i])
        if isempty(isozymes)
            push!(columns, _gecko_column(i, 0, 0, 0, lbs[i], ubs[i], [], []))
            continue
        end

        # if the reaction has multiple isozymes, it needs extra coupling to
        # ensure that the total rate of the reaction doesn't exceed the
        # "global" limit
        if relaxed_arm_reaction_bounds
            reaction_coupling_row =
                length(isozymes) > 1 ? begin
                    push!(coupling_row_reaction, i)
                    length(coupling_row_reaction)
                end : 0
        end

        # loop over both directions for all isozymes
        for (lb, ub, kcatf, dir) in [
            (-ubs[i], -lbs[i], i -> i.kcat_reverse, -1),
            (lbs[i], ubs[i], i -> i.kcat_forward, 1),
        ]
            if !relaxed_arm_reaction_bounds
                # In this case, the coefficients in the coupling matrix will be
                # the same as in the combined case, only categorized in
                # separate rows for negative and positive ones. Surprisingly,
                # we do not need to explicitly remember the bounds, because the
                # ones taken from the original model are perfectly okay -- the
                # "reverse" direction is unreachable because of individual
                # bounds on split reactions, and the "forward" direction is
                # properly negated in the reverse case to work nicely with the
                # global lower bound.
                reaction_coupling_row =
                    ub > 0 && length(isozymes) > 1 ? begin
                        push!(coupling_row_reaction, i)
                        length(coupling_row_reaction)
                    end : 0
            end

            # all isozymes in this direction
            for (iidx, isozyme) in enumerate(isozymes)
                kcat = kcatf(isozyme)
                if ub > 0 && kcat > _constants.tolerance
                    # prepare the coupling with gene product molar
                    gene_product_coupling = collect(
                        begin
                            gidx = gene_name_lookup[gene]
                            row_idx = if haskey(gene_row_lookup, gidx)
                                gene_row_lookup[gidx]
                            else
                                push!(coupling_row_gene_product, gidx)
                                gene_row_lookup[gidx] =
                                    length(coupling_row_gene_product)
                            end
                            (row_idx, stoich / kcat)
                        end for (gene, stoich) in isozyme.gene_product_count if
                        haskey(gene_name_lookup, gene)
                    )

                    # prepare the coupling with the mass groups
                    gp_groups = gmg_.(keys(isozyme.gene_product_count))
                    gp_mass = gpmm_.(keys(isozyme.gene_product_count))
                    groups = unique(filter(!isnothing, gp_groups))
                    group_idx = Dict(groups .=> 1:length(groups))
                    vals = [0.0 for _ in groups]

                    for (gpg, mass) in zip(gp_groups, gp_mass)
                        if !isnothing(gpg)
                            vals[group_idx[gpg]] += mass / kcat
                        end
                    end

                    mass_group_coupling = collect(
                        isnothing(group) ? 0 :
                        begin
                            if !haskey(mass_group_lookup, group)
                                push!(coupling_row_mass_group, group)
                                mass_group_lookup[group] =
                                    length(coupling_row_mass_group)
                            end
                            (mass_group_lookup[group], val)
                        end for (group, val) in zip(groups, vals)
                    )

                    # make a new column
                    push!(
                        columns,
                        _gecko_column(
                            i,
                            iidx,
                            dir,
                            reaction_coupling_row,
                            max(lb, 0),
                            ub,
                            gene_product_coupling,
                            mass_group_coupling,
                        ),
                    )
                end
            end
        end
    end

    GeckoModel(
        columns,
        coupling_row_reaction,
        collect(zip(coupling_row_gene_product, gpb_.(gids[coupling_row_gene_product]))),
        collect(zip(coupling_row_mass_group, gmgb_.(coupling_row_mass_group))),
        model,
    )
end
