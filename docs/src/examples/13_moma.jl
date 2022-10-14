# # Minimization of metabolic adjustment (MOMA)

# MOMA allows you to find a feasible solution of the model that is closest (in
# an Euclidean metric) to a reference solution. Often this gives a realistic
# estimate of the organism behavior that has undergone a radical change (such
# as a gene knockout) that prevents it from metabolizing optimally, but the
# rest of the metabolism has not yet adjusted to compensate for the change.

# The original description of MOMA is by [Segre, Vitkup, and Church, "Analysis
# of optimality in natural and perturbed metabolic networks", Proceedings of the
# National Academy of Sciences, 2002](https://doi.org/10.1073/pnas.232349399).

# As always, let's start with downloading a model.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA

model = load_model(StandardModel, "e_coli_core.xml")

# MOMA analysis requires solution of a quadratic model, we will thus use Clarabel as the main optimizer.

using Clarabel

# We will need a reference solution, which represents the original state of the
# organism before the change.
reference_flux =
    flux_balance_analysis_dict(model, Clarabel.Optimizer; modifications = [silence])

# As the change here, we manually knock out CYTBD reaction:
changed_model = change_bound(model, "R_CYTBD", lower = 0.0, upper = 0.0);

# Now, let's find a flux that minimizes the organism's metabolic adjustment for
# this model:
flux_summary(
    minimize_metabolic_adjustment_analysis_dict(
        changed_model,
        reference_flux,
        Clarabel.Optimizer;
        modifications = [silence],
    ),
)

# For illustration, you can compare the result to the flux that is found by
# simple optimization:

flux_summary(
    flux_balance_analysis_dict(
        changed_model,
        Clarabel.Optimizer;
        modifications = [silence],
    ),
)
