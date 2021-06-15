# The "full" analysis:
# all regimes (min, median, max, and Miura)
# no duplicate proviruses in the integration dates
# special handling for p3

load("full_analysis_ode.RData")

source("analysis_helpers.r")

# Read in the sample time data.
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
    # For this analysis, we remove duplicate sequences.
    subject.data <- subject.data[is.na(subject.data$duplicate),]

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

save(
    all.subjects,
    bin.size,
    possible.half.lives,
    all.log.likelihoods,
    mles,
    bayes.factors,
    file="full_analysis.RData"
)
