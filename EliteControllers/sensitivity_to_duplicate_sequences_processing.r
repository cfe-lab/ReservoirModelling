# This analysis uses:
# - standard integration dates
# - duplicate proviruses *included*
# - Miura regime only
# - special handling for P3's viral load data

# This is commented out because we'll typically already have loaded the data.
# load("full_analysis_ode.RData")

source("analysis_helpers.r")

# We only consider the Miura regime.
acute.phase[["min"]] <- NULL
acute.phase[["median"]] <- NULL
acute.phase[["max"]] <- NULL

regimes <- "Miura"
subjects <- c("p1", "p2", "p3", "p4")

for (subject in subjects) {
    vl.data[[subject]][["min"]] <- NULL
    vl.data[[subject]][["median"]] <- NULL
    vl.data[[subject]][["max"]] <- NULL
}

# The ODEs are unaffected by the different integration date data, but we
# can eliminate everything for regimes except for the Miura regime.
for (bin.size.label in c("bin.30", "bin.365")) {
    for (subject in subjects) {
        ode.solutions[[bin.size.label]][[subject]][["min"]] <- NULL
        ode.solutions[[bin.size.label]][[subject]][["median"]] <- NULL
        ode.solutions[[bin.size.label]][[subject]][["max"]] <- NULL
    }
}

integration.data <- prepare.integration.data(
    subjects,
    vl.data,
    remove.duplicates=FALSE,
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
    integration.data,
    ode.solutions$bin.30,
    bin.size,
    possible.half.lives
)

bayes.factors <- compute.bayes.factors(
    integration.data,
    ode.solutions$bin.30,
    all.log.likelihoods
)

# Here, we save the whole kaboodle because we pared down the ODEs too.
save.image("sensitivity_to_duplicate_sequences.RData")
