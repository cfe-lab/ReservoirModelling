# May 3, 2021: Zabrina spotted an issue with the non-Ragon controllers analysis
# that revealed that both that analysis and this one had incorrect values for the VL:
# they were specified in copies/mL and the ODE expects copies/uL.  This script
# redoes the analyses with the viral load corrected.

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
    sample.times.base.path <- "../../data"
    if (subject == "p3") {
        sample.times.base.path <- paste(sample.times.base.path, "EliteControllers_2021_03_01", sep="/")
    }
    # The sample times CSV file has 4 useful columns, toss the rest.
    subject.data <- read.csv(paste(sample.times.base.path, subject, "sample_times.csv", sep="/"))
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

save.image("full_analysis.RData")
