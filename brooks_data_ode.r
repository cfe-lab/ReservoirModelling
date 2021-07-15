# Solve the ODEs numerically (this is the slow part).
# Note that this is completely independent of the integration data.

source("ode_based_helpers.r")
source("reservoir_helpers.r")

vl.data <- load.vl.data()
ode.results <- perform.odes(vl.data)

save(
    ode.results,
    file="brooks_data_ode.RData"
)
