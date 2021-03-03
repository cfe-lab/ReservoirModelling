# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.

# This is commented out as we will often have already loaded it when we run this script.
# load("full_analysis_2021_03_03.RData")


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
        breakpoints = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=365)
        if (!(max(curr.data$days.before.art) %in% breakpoints)) {
            breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + 365)
        }
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

        # Make some plots of the regressions.
        fit <- predict(decay.rate.regression, se.fit=TRUE)

        predictions <- exp(fit$fit)
        upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
        lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)

        max.y <- max(c(actual.freqs$counts, upper.bounds))

        pdf(paste("model_free_decay_rate_", subject, "_", regime, "_known_vl.pdf", sep=""))
        par(mar=c(6.5, 6.5, 2, 2) + 0.1)
        plot(
            c(0, length(actual.freqs$counts)),
            c(0, max.y),
            xlab=NA,
            ylab=NA,
            type="n",
            xaxt="n",
            cex.lab=3,
            cex.axis=2
        )

        # Some defaults, and then some customization for p3.
        x.label <- "Year prior to ART initiation"
        x.label.cex <- 3
        if (subject == "p3") {
            x.label <- "Year prior to last viremic episode"
            x.label.cex <- 2.25
        }

        title(
            xlab=x.label,
            line=3.5,
            cex.lab=x.label.cex
        )

        title(
            ylab="Provirus proportion",
            line=3.5,
            cex.lab=3
        )

        # Plot the actual observed counts.
        lines(
            actual.freqs$counts,
            type="h",
            lwd=20,
            col="blue",
            lend="butt"
        )

        # Overlay the GLM fitted value and error bars.
        points(predictions)
        arrows(
            x0=1:length(actual.freqs$counts),
            y0=lower.bounds,
            y1=upper.bounds,
            code=3,
            angle=90
        )
        dev.off()
    }
}
