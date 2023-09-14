# # Hit and run sampling

# Sampling the feasible space of the model allows you to gain a realistic
# insight into the distribution of the flow and its probabilistic nature, often
# better describing the variance and correlations of the possible fluxes better
# (but more approximately) than e.g. [flux variability analysis](06_fva.md).

# COBREXA supports a variant of hit-and-run sampling adjusted to the
# complexities of metabolic models; in particular, it implements a version
# where the next sample (and indirectly, the next run direction) is generated
# as an affine combination of the samples in the current sample set. This gives
# a result similar to the artificially centered hit-and-run sampling, with a
# slightly (unsubstantially) different biases in the generation of the next
# hits.

# As always, we start by loading everything that is needed:

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, GLPK

model = load_model("e_coli_core.xml")

# The sampling procedure requires a set of "seed" points that will form the
# basis for the first iteration of runs. This is commonly called a warm-up. You
# can generate these points manually with any method of choice, or you can use
# the COBREXA implementation of warmup that uses the extreme edges of the
# polytope, similarly to the way used by flux variability analysis:

warmup_points = warmup_from_variability(model, GLPK.Optimizer)

# This generates a matrix of fluxes, where each column is one sample point, and
# rows correspond to the reactions in the model.

# With this, starting the sampling procedure is straightforward:

samples = affine_hit_and_run(model, warmup_points)

# As seen above, the result again contains 1 sample per column, with reactions
# in the same order with rows as before. To get a different sample, you can
# tune the parameters the function. `sample_iters` allows you to specify the
# iterations at which you want the current sample to be collected and reported
# -- by default, that is done 5 times on each 100th iteration. In the example
# below, we catch the samples in 10 iterations right after the 200th iteration
# passed.  Similarly, to avoid possible degeneracy, you can choose to run more
# hit-and-run batch chains than 1, using the `chains` parameters. The total
# number of points collected is the number of points in warmup times the number
# of sample-collecting iterations times the number of chains.

samples = affine_hit_and_run(model, warmup_points, sample_iters = 201:210, chains = 2)

#md # !!! tip "Parallelization"
#md #     Both procedures used for sampling in this example
#md #     ([`warmup_from_variability`](@ref), [`affine_hit_and_run`](@ref)) can be
#md #     effectively parallelized by adding `workers=` parameter, as summarized
#md #     in [the documentation](../distributed/1_functions.md). Due to the nature of the algorithm, parallelization
#md #     of the sampling requires at least 1 chain per worker.

# ## Visualizing the samples

# Samples can be displayed very efficiently in a scatterplot or a density plot,
# which naturally show correlations and distributions of the fluxes:

using CairoMakie

o2, co2 = indexin(["R_EX_o2_e", "R_EX_co2_e"], variable_ids(model))

scatter(
    samples[o2, :],
    samples[co2, :];
    axis = (; xlabel = "Oxygen exchange", ylabel = "Carbon dioxide exchange"),
)
