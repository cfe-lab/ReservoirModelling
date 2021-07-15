# Estimate the decay rates via MLE *with duplicate proviruses included*.

source("reservoir_helpers.r")
source("ode_based_helpers.r")

# Note: we only use this to do a sanity check on our data.
vl.data <- load.vl.data()

# This might already have been done.
# load("brooks_data_ode.RData")

without.duplicates <- load.integration.data()
with.duplicates <- load.integration.data.with.duplicates()

# Sanity check: check that for all of the PIDs represented in the 
# "with duplicates" data, the contents of with.duplicates is equal
# to that in without.duplicates.
for (pid in unique(with.duplicates$pid)) {
    curr.without <- without.duplicates[without.duplicates$pid == pid, c(1:6, 9:11)]
    curr.with <- with.duplicates[with.duplicates$pid == pid,]

    equality.check <- curr.with[1:nrow(curr.without),] == curr.without  # this is a matrix
    if (!all(equality.check)) {
        cat(
            "Warning: there appears to be a discrepancy in the \"with duplicates\" ",
            "data for PID ",
            pid,
            "\n",
            sep=""
        )
    }
}

integration.data <- load.combined.integration.data(
    without.duplicates,
    with.duplicates
)

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
