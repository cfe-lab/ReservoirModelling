# Estimate the decay rates via MLE.

source("reservoir_helpers.r")
source("ode_based_helpers.r")

# Note: we only use this to do a sanity check on our data.
vl.data <- load.vl.data()

# This might already have been done.
# load("brooks_data_ode.RData")

integration.data <- load.integration.data()

date.sanity.check(vl.data, integration.data)

# These are default values:
bin.size <- 30
possible.half.lives <- (1:200) * bin.size

# Find the decay rate that maximizes the likelihood for each individual.
estimation.results <- ode.based.likelihood.estimates(
    integration.data,
    ode.results$bin.30,
    bin.size=bin.size,
    possible.half.lives=possible.half.lives
)

all.log.likelihoods <- estimation.results$all.log.likelihoods
mles <- estimation.results$mles
