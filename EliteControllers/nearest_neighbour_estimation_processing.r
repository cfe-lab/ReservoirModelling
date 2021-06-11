# This is a similar analysis to the main analysis, but using integration
# dates estimated using the nearest-neighbour method.

source("analysis_helpers.r")

# This is commented out because we will typically have already loaded this,
# or have just re-run the analysis producing it.
# load("full_analysis.RData")

# Read in the sample time data.
for (subject in subjects) {
    subject.data <- read.csv(
        paste(
            "../../data/NearestNeighbourEstimation_2021_06_08", 
            paste0(
                subject,
                "_integration.csv"
            ),
            sep="/"
        )
    )
    # These files have two useful columns (the first and second):
    # PBMC ID
    # NN Year (which is actually a date in YYYY-MM-DD format)
    subject.data <- subject.data[, c(1, 2)]
    names(subject.data) <- c(
        "id",
        "integration.date.est"
    )
    subject.data$integration.date.est <- strptime(subject.data$integration.date.est, "%Y-%m-%d")

    # No special handling is needed for P3's integration dates, as with this method
    # they're all necessarily before cART initiation.
    subject.data$days.before.art <- compute.days.before.art(
        subject.data$integration.date.est,
        all.subjects[[subject]]$art.initiation,
        all.subjects[[subject]]$infection.date
    )

    # Overwrite the integration data that we inherit from the full analysis.
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

save.image("nn_integration_dates_analysis.RData")
