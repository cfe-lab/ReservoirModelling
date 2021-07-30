# This analysis is:
# p3 only
# no special handling on P3's integration dates
# including duplicate proviruses

# This will typically already be loaded.
# load("p3_no_special_handling_ode.RData")

source("analysis_helpers.r")

# Unpack the results of the loaded preamble.
vl.data <- results.p3.no.special.handling$read.vl$vl.data
ode.solutions <- results.p3.no.special.handling$ode.solutions

integration.data <- prepare.integration.data(
    "p3",
    vl.data,
    remove.duplicates=FALSE,
    p3.boundaries=NULL
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
