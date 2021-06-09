# May 4, 2021: redoing the analysis of March 3 but with the viral load corrected.

source("analysis_helpers.r")

# These constants are taken from the literature.
acute.phase <- list()
acute.phase[["Miura"]] <- list(
    days.to.peak=31,
    peak.vl=53300,
    days.to.undetectable=41
)
regimes <- names(acute.phase)
subjects <- c("p1", "p2", "p3", "p4")

# Read in the viral load data.
read.vl <- prepare.vl.data(subjects, p3.special.handling=TRUE)
all.subjects <- read.vl$all.subjects
# Some special dates that we will have to remember for p3.
p3.art.initiation <- read.vl$p3.art.initiation
p3.blip <- read.vl$p3.blip

# Graft on the acute phase.
grafted.vls <- graft.acute.phase(all.subjects, acute.phase, acute.setpoint=1000)

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


# Solve the ODEs numerically (this is the slow part).
ode.solutions <- solve.odes(all.subjects, grafted.vls)
ode.solutions.bin.30 <- ode.solutions$bin.30
ode.solutions.bin.365 <- ode.solutions$bin.365


# Find the decay rate that maximizes the likelihood for each individual.
bin.size <- 30
possible.half.lives <- (1:1800) * bin.size  # months, roughly
all.log.likelihoods <- compute.lls(
    all.subjects,
    ode.solutions$bin.30,
    bin.size=bin.size,
    possible.half.lives=possible.half.lives
)

save.image("sensitivity_to_duplicate_sequences.RData")
