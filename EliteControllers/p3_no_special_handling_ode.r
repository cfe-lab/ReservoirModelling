# This is preamble for analysis that cover:
# p3 only
# no special handling of p3's VL data

# We start from scratch so that we can redo the ODE with the un-adjusted VLs.

source("analysis_helpers.r")

# Read in the viral load data.
read.vl <- prepare.vl.data("p3", p3.special.handling=FALSE)
vl.data <- read.vl$vl.data

# Graft on the acute phase.
grafted.vls <- graft.acute.phase(vl.data, acute.phase, acute.setpoint=1000)

# Solve the ODEs numerically (this is the slow part).
ode.solutions <- solve.odes(vl.data, grafted.vls)

results.p3.no.special.handling <- list(
    read.vl=read.vl,
    grafted.vls=grafted.vls,
    ode.solutions=ode.solutions
)

save(results.p3.no.special.handling, file="p3_no_special_handling_ode.RData")
