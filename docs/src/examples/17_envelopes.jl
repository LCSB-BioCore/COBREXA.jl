
# # Production envelopes

# Production envelopes show what a model is capable of doing on a wide range of
# parameters. Usually, you choose a regular grid of a small dimension in the
# parameter space, and get an information about how well the model runs at each
# point.

# As usual, we start by loading everything:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model("e_coli_core.xml")

# The envelope analysis "structure" is similar to what you can obtain using
# [`screen`](@ref); but COBREXA provides a special functions that run the
# process in a very optimized manner. For envelopes, there is
# [`envelope_lattice`](@ref) that generates the rectangular lattice of points
# for a given model and reactions, and [`objective_envelope`](@ref) that
# computes the output (usually as the objective value) of the model at the
# lattice points. You do not need to call [`envelope_lattice`](@ref) directly
# because it is taken as a "default" way to create the lattice by
# [`objective_envelope`](@ref).

# In short, we can compute the envelope of a single reaction in the *E. coli*
# model as follows:

envelope = objective_envelope(
    model,
    ["R_EX_o2_e"],
    GLPK.Optimizer,
    lattice_args = (ranges = [(-50, 0)],),
)

# (The named tuple given in `lattice_args` argument is passed to the internal
# call of [`envelope_lattice`](@ref), giving you an easy way to customize its
# behavior.)

# The result has 2 fields which can be used to easily plot the envelope. We
# also need to "fix" the missing values (represented as `nothing`) where the
# model failed to solve -- we will simply omit them here).

using CairoMakie

valid = .!(isnothing.(envelope.values))
lines(envelope.lattice[1][valid], float.(envelope.values[valid]))

# Additional resolution can be obtained either by supplying your own larger
# lattice, or simply forwarding the `samples` argument to the internal call of
# [`envelope_lattice`](@ref):

envelope = objective_envelope(
    model,
    ["R_EX_co2_e"],
    GLPK.Optimizer,
    lattice_args = (samples = 1000, ranges = [(-50, 100)]),
)

valid = .!(isnothing.(envelope.values))
lines(envelope.lattice[1][valid], float.(envelope.values[valid]))

# ## Multi-dimensional envelopes

# The lattice functions generalize well to more dimensions; you can easily
# explore the production of the model depending on the relative fluxes of 2
# reactions:

envelope = objective_envelope(
    model,
    ["R_EX_o2_e", "R_EX_co2_e"],
    GLPK.Optimizer,
    lattice_args = (samples = 100, ranges = [(-60, 0), (-15, 60)]),
)

heatmap(
    envelope.lattice[1],
    envelope.lattice[2],
    [isnothing(x) ? 0 : x for x in envelope.values],
    axis = (; xlabel = "Oxygen exchange", ylabel = "Carbon dioxide exchange"),
)
