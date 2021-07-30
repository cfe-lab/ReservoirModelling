# The "full" analysis:
# all regimes (min, median, max, and Miura)
# no duplicate proviruses in the integration dates
# special handling for p3

# We comment this out as this is typically already loaded.
# load("full_analysis_ode.RData")

source("analysis_helpers.r")

integration.data <- prepare.integration.data(
    subjects,
    vl.data,
    remove.duplicates=TRUE,
    p3.boundaries=c(p3.art.initiation, p3.blip)
)


# Find the decay rate that maximizes the likelihood for each individual.
bin.size <- 30
possible.half.lives <- (1:1800) * bin.size  # months, roughly
all.log.likelihoods <- compute.lls(
    integration.data,
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
    integration.data,
    ode.solutions$bin.30,
    all.log.likelihoods
)

save(
    integration.data,
    bin.size,
    possible.half.lives,
    all.log.likelihoods,
    mles,
    bayes.factors,
    file="full_analysis.RData"
)
