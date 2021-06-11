# May 4, 2021: redoing the analysis of March 3 but with the viral load corrected.

# This is commented out because we'll typically already have loaded the data.
# load("full_analysis.RData")

source("analysis_helpers.r")

# We only consider the Miura regime.
acute.phase[["min"]] <- NULL
acute.phase[["median"]] <- NULL
acute.phase[["max"]] <- NULL

regimes <- "Miura"
subjects <- c("p1", "p2", "p3", "p4")

for (subject in subjects) {
    all.subjects[[subject]][["min"]] <- NULL
    all.subjects[[subject]][["median"]] <- NULL
    all.subjects[[subject]][["max"]] <- NULL
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


# Read in the sample time data, this time retaining the duplicates.
for (subject in subjects) {
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

    if (subject != "p3") {
        subject.data$days.before.art <- compute.days.before.art(
            subject.data$integration.date.est,
            all.subjects[[subject]]$art.initiation,
            all.subjects[[subject]]$infection.date
        )
    } else {
        subject.data$days.before.art <- compute.days.before.art(
            subject.data$integration.date.est,
            all.subjects[[subject]]$art.initiation,
            all.subjects[[subject]]$infection.date,
            boundaries=c(p3.art.initiation, p3.blip)
        )
    }
    all.subjects[[subject]][["integration"]] <- subject.data
}


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

save.image("sensitivity_to_duplicate_sequences.RData")
