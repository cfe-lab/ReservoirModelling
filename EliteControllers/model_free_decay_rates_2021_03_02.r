# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.

# This is commented out as we will often have already loaded it when we run this script.
# load("elite_controllers_known_vl_2021_03_01.RData")


# Plots of reservoir composition.
all.regressions <- list()
for (subject in subjects) {
    all.regressions[[subject]] <- list()
    curr.data <- all.subjects[[subject]]$integration

    days.pre.therapy <- 
        as.numeric(
            all.subjects[[subject]]$art.initiation - all.subjects[[subject]]$infection.date,
            units="days"
        )

    for (regime in regimes) {
        all.regressions[[subject]][[regime]] <- list()
        breakpoints = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=365)
        if (!(max(curr.data$days.before.art) %in% breakpoints)) {
            breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + 365)
        }
        actual.freqs <- hist(
            curr.data$days.before.art,
            breaks=breakpoints,
            plot=FALSE
        )

        counts.by.year <- rev(actual.freqs$counts)
        decay.rate.regression <- glm(
            actual.freqs$counts ~ 1:length(counts.by.year),
            family=poisson
        )
        all.regressions[[subject]][[regime]] <- decay.rate.regression
    }
}
