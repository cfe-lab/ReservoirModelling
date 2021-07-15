# Half-life plots with estimates drawn from the "with-model" calculations.

source("ode_based_helpers.r")

# This might already have been done.
# load("brooks_data_ode.RData")

integration.data <- load.integration.data()

# Find the decay rate that maximizes the likelihood for each individual.
estimation.results <- ode.based.likelihood.estimates(
    integration.data,
    ode.results$bin.30
)
mles <- estimation.results$mles

# Our estimates are in the data frame `mles`, which we need to refactor
# to plug into our plotting methods: these estimates are in days and we
# need to convert them to years.
mles$half.life <- mles$mle / 365
mles$lower.bound <- mles$lower.bound / 365
mles$upper.bound <- mles$upper.bound / 365

source("half_life_plot_helpers.r")

pdf("ode_based_known_vl_half_lives.pdf")
plot.all.individuals.half.lives(mles[mles$vl.info == "known",])
dev.off()

pdf("ode_based_known_vl_half_lives_two_sampling_points.pdf")
plot.multiple.timepoints.half.lives(mles[mles$vl.info == "known",])
dev.off()

pdf("ode_based_typical_vl_half_lives.pdf")
plot.all.individuals.half.lives(mles[mles$vl.info == "typical",])
dev.off()

pdf("ode_based_typical_vl_half_lives_two_sampling_points.pdf")
plot.multiple.timepoints.half.lives(mles[mles$vl.info == "typical",])
dev.off()
