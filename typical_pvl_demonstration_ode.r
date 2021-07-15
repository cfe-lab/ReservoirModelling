source("reservoir_helpers.r")

cat("Calculating undecayed distribution of reservoir....\n")
run.time <- system.time(
    ode.bin.365 <- undecayed.reservoir.distribution(
        as.numeric(10 * 365, units="days"),
        bin.size=365
    )
)
cat("Time elapsed:\n")
print(run.time)

save.image("typical_pvl_demonstration_ode.RData")
