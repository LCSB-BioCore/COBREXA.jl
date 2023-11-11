import JSONFBCModels as M
import AbstractFBCModels as A
import COBREXA as X
import ConstraintTrees as C
import JuMP as J
import Gurobi as G

model = A.load(M.JSONFBCModel, "e_coli_core.json")

reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}()
for rid in A.reactions(model)
    haskey(ecoli_core_reaction_kcats, rid) || continue
    for (i, grr) in enumerate(A.reaction_gene_association_dnf(model, rid))
        d = get!(reaction_isozymes, rid, Dict{String,X.Isozyme}())
        d["isozyme_"*string(i)] = X.Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))),
            kcat_forward = ecoli_core_reaction_kcats[rid][i][1],
            kcat_backward = ecoli_core_reaction_kcats[rid][i][2],
        )
    end
end

m = X.metabolic_model(model)

# create directional fluxes
m +=
    :fluxes_forward^X.fluxes_in_direction(m.fluxes, :forward) +
    :fluxes_backward^X.fluxes_in_direction(m.fluxes, :backward)

# link directional fluxes to original fluxes
m *=
    :link_flux_directions^X.sign_split_constraints(
        positive = m.fluxes_forward,
        negative = m.fluxes_backward,
        signed = m.fluxes,
    )

# create fluxes for each isozyme
for (rid, _) in m.fluxes_forward
    if haskey(reaction_isozymes, string(rid))
        m += :fluxes_isozymes_forward^rid^X.isozyme_variables(string(rid), reaction_isozymes)
    end
end
for (rid, _) in m.fluxes_backward
    if haskey(reaction_isozymes, string(rid))
        m += :fluxes_isozymes_backward^rid^X.isozyme_variables(string(rid), reaction_isozymes)
    end
end

# link isozyme fluxes to directional fluxes
m *=
    :link_isozyme_fluxes_forward^X.link_isozymes(
        m.fluxes_forward,
        m.fluxes_isozymes_forward,
    )
m *=
    :link_isozyme_fluxes_backward^X.link_isozymes(
        m.fluxes_backward,
        m.fluxes_isozymes_backward,
    )

# create enzyme variables
m += :enzymes^X.enzyme_variables(model)

# add enzyme mass balances
m *=
    :enzyme_stoichiometry^X.enzyme_stoichiometry(
        m.enzymes,
        m.fluxes_isozymes_forward,
        m.fluxes_isozymes_backward,
        reaction_isozymes,
    )

# add capacity limitations
m *=
    :capacity_limitation^X.enzyme_capacity(
        m.enzymes,
        enzyme_molar_mass,
        [:G_b0351, :G_b1241, :G_b0726, :G_b0727],
        10.0,
    )
