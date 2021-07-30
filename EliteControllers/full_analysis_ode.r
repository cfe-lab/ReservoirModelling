# The idea behind this script is to put all the slow stuff right at the beginning,
# so that all subsequent analyses can simply load the result.

source("analysis_helpers.r")

# Read in the viral load data.
read.vl <- prepare.vl.data(subjects, p3.special.handling=TRUE)
vl.data <- read.vl$vl.data
# Some special dates that we will have to remember for p3.
p3.art.initiation <- read.vl$p3.art.initiation
p3.blip <- read.vl$p3.blip

# Graft on the acute phase.
grafted.vls <- graft.acute.phase(vl.data, acute.phase, acute.setpoint=1000)

# Solve the ODEs numerically (this is the slow part).
ode.solutions <- solve.odes(vl.data, grafted.vls)

save.image("full_analysis_ode.RData")
