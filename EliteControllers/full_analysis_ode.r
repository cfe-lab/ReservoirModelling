# The idea behind this script is to put all the slow stuff right at the beginning,
# so that all subsequent analyses can simply load the result.

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

# Solve the ODEs numerically (this is the slow part).
ode.solutions <- solve.odes(all.subjects, grafted.vls)

save.image("full_analysis_ode.RData")
