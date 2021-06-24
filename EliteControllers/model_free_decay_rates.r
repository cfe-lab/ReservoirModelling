# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.

# This is commented out as we will often have already loaded it when we run this script.
# load("full_analysis_ode.RData")
# load("full_analysis.RData")

source("plot_helpers.r")

all.regressions <- list()
for (subject in subjects) {
    all.regressions[[subject]] <- list()
    curr.data <- all.subjects[[subject]]$integration

    days.pre.therapy <- 
        as.numeric(
            all.subjects[[subject]]$art.initiation - all.subjects[[subject]]$infection.date,
            units="days"
        )

    for (regime in "Miura") {
        all.regressions[[subject]][[regime]] <- list()

        breakpoints <- compute.bin.breakpoints(
            curr.data$days.before.art,
            days.pre.therapy,
            bin.size=365
        )
        actual.freqs <- hist(
            curr.data$days.before.art,
            breaks=breakpoints,
            plot=FALSE
        )
        regression.frame <- data.frame(
            x=1:length(actual.freqs$counts),
            y=actual.freqs$counts
        )
        decay.rate.regression <- glm(
            y ~ x,
            family=poisson,
            data=regression.frame
        )
        all.regressions[[subject]][[regime]] <- decay.rate.regression

        # As per this link:
        # http://biometry.github.io/APES/Stats/stats23-GeneralizedLinearModels-GLM.html
        # there is *some* evidence of overdispersion in this fit, but not as decisive
        # as in the example.

        x.label="Year prior to ART initiation"
        if (subject == "p3") {
            x.label="Year prior to last viremic episode"
        }
        cairo_pdf(paste("model_free_decay_rate_", subject, ".pdf", sep=""))
        model.free.plot(
            decay.rate.regression,
            actual.freqs$counts,
            x.label=x.label
        )
        dev.off()
    }
}
