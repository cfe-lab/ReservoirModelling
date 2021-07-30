# This is a similar analysis to the main analysis, but using integration
# dates estimated using the nearest-neighbour method.

source("analysis_helpers.r")

# This is commented out because we will typically have already loaded this,
# or have just re-run the analysis producing it.
# load("full_analysis_ode.RData")

# Read in the sample time data.
nn.integration.data <- prepare.nearest.neighbour.integration.data(
    subjects,
    vl.data
)

# Find the decay rate that maximizes the likelihood for each individual.
bin.size <- 30
possible.half.lives <- (1:1800) * bin.size  # months, roughly
all.log.likelihoods <- compute.lls(
    nn.integration.data,
    ode.solutions$bin.30,
    bin.size=bin.size,
    possible.half.lives=possible.half.lives
)

mles <- compute.mles(
    all.log.likelihoods,
    ode.solutions$bin.30,
    bin.size,
    possible.half.lives
)

bayes.factors <- compute.bayes.factors(
    nn.integration.data,
    ode.solutions$bin.30,
    all.log.likelihoods
)

save(
    nn.integration.data,
    bin.size,
    possible.half.lives,
    all.log.likelihoods,
    mles,
    bayes.factors,
    file="nearest_neighbour_estimation_analysis.RData"
)
