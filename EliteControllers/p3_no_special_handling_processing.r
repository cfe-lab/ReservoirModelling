# This analysis is:
# p3 only
# no special handling on P3's integration dates
# no duplicate proviruses
# We have to start from scratch so that we can redo the ODE with the un-adjusted 

source("analysis_helpers.r")

# These constants are taken from the literature.
acute.phase <- list()
acute.phase[["min"]] <- list(
    days.to.peak=14,
    peak.vl=1954,
    days.to.undetectable=6 * 7
)
acute.phase[["median"]] <- list(
    days.to.peak=23,
    peak.vl=12902,
    days.to.undetectable=22 * 7
)
acute.phase[["max"]] <- list(
    days.to.peak=30,
    peak.vl=71550,
    days.to.undetectable=43 * 7
)
acute.phase[["Miura"]] <- list(
    days.to.peak=31,
    peak.vl=53300,
    days.to.undetectable=41
)
regimes <- names(acute.phase)
subjects <- "p3"
subject <- "p3"

# Read in the viral load data.
read.vl <- prepare.vl.data(subjects, p3.special.handling=FALSE)
all.subjects <- read.vl$all.subjects

# Graft on the acute phase.
grafted.vls <- graft.acute.phase(all.subjects, acute.phase, acute.setpoint=1000)

# Solve the ODEs numerically (this is the slow part).
ode.solutions <- solve.odes(all.subjects, grafted.vls)

# Eliminate several columns we don't use.
subject.data <- read.csv(
    paste(
        "../../data/IntegrationData_2021_06_10",
        paste0(
            subject,
            "_integration.csv"
        ),
        sep="/"
    )
)
subject.data <- subject.data[, c(1, 4, 5, 6, 7)]
names(subject.data) <- c(
    "id",
    "integration.date.est",
    "integration.date.lower",
    "integration.date.upper",
    "duplicate"
)
for (col.idx in 2:4) {
    subject.data[[col.idx]] <- strptime(subject.data[[col.idx]], "%Y-%m-%d")
}

# No special handling for the integration dates.
subject.data$days.before.art <- compute.days.before.art(
    subject.data$integration.date.est,
    all.subjects[["p3"]]$art.initiation,
    all.subjects[["p3"]]$infection.date
)

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

# Save the whole kaboodle as we pared down the ODEs too.
save.image("p3_no_special_handling_analysis.RData")
