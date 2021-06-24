# This analysis is:
# p3 only
# no special handling on P3's integration dates
# no duplicate proviruses

source("analysis_helpers.r")

# This will typically already be loaded.
# load("p3_no_special_handling_including_duplicates_analysis.RData")

# We inherit nearly everything from the above.  Here, we just need
# to remove duplicate sequences.
subject.data <- subject.data[is.na(subject.data$duplicate),]
all.subjects[["p3"]][["integration"]] <- subject.data


# Find the decay rate that maximizes the likelihood for each individual.
bin.size <- 30
possible.half.lives <- (1:1800) * bin.size  # months, roughly
all.log.likelihoods <- compute.lls(
    all.subjects,
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
    all.subjects,
    ode.solutions$bin.30,
    all.log.likelihoods
)
