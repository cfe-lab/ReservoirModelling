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

# Read in the sample time data.
for (subject in subjects) {
    sample.times.path <- ""
    # This is where the proviral integration data *with duplicates* lives.
    if (subject != "p3") {
        sample.times.path <- paste(
            "../../data/EliteControllers_2020_11_30",
            subject,
            "sample_times.csv",
            sep="/"
        )
    } else {
        sample.times.path <- "../../data/P3_with_duplicates_no_truncation.csv"
    }
    subject.data <- read.csv(sample.times.path)

    # The sample times CSV file has 4 useful columns, toss the rest.
    subject.data <- subject.data[, c(1, 4, 5, 6)]
    names(subject.data) <- c(
        "id",
        "integration.date.est",
        "integration.date.lower",
        "integration.date.upper"
    )
    for (col.idx in 2:4) {
        subject.data[[col.idx]] <- strptime(subject.data[[col.idx]], "%Y-%m-%d")
    }

    if (subject != "p3") {
        subject.data$days.before.art <- compute.days.before.art(
            subject.data$integration.date.est,
            all.subjects[[subject]]$art.initiation
        )
    } else {
        subject.data$days.before.art <- compute.days.before.art(
            subject.data$integration.date.est,
            all.subjects[[subject]]$art.initiation,
            boundaries=c(p3.art.initiation, p3.blip)
        )
    }
    all.subjects[[subject]][["integration"]] <- subject.data
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

save.image("sensitivity_to_duplicate_sequences.RData")
