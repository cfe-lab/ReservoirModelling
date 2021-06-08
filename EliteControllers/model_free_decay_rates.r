# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.

# This is commented out as we will often have already loaded it when we run this script.
# load("full_analysis.RData")


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

        # As per this link:
        # http://biometry.github.io/APES/Stats/stats23-GeneralizedLinearModels-GLM.html
        # there is *some* evidence of overdispersion in this fit, but not as decisive
        # as in the example.

        # Make some plots of the regressions.
        # The error bars were computed with the guidance of this helpful blog post: 
        # https://www.r-bloggers.com/2018/12/confidence-intervals-for-glms/
        x.fit.values <- data.frame(x=seq(0.5, length(actual.freqs$counts) + 0.5, by=0.1))
        fit <- predict(decay.rate.regression, newdata=x.fit.values, se.fit=TRUE)

        predictions <- exp(fit$fit)
        upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
        lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)

        total.count <- sum(actual.freqs$counts)
        # max.y <- max(c(actual.freqs$counts, upper.bounds)) / total.count
        # # Some customization for p2.
        # if (subject == "p2") {
        #     max.y <- max.y + 0.075
        # }
        # max.y <- min(1, max.y)
        max.y <- 1

        cairo_pdf(paste("model_free_decay_rate_", subject, ".pdf", sep=""))
        par(mar=c(6.5, 6.5, 2, 2) + 0.1)
        plot(
            c(0, length(actual.freqs$counts)),
            c(0, max.y),
            xlab=NA,
            ylab=NA,
            type="n",
            xaxt="n",
            yaxt="n",
            cex.lab=3,
            cex.axis=2
        )

        # Some defaults, and then some customization for p3.
        x.label <- "Year prior to ART initiation"
        x.label.cex <- 2.25
        if (subject == "p3") {
            x.label <- "Year prior to last viremic episode"
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

        axis(
            1,  # this is the x axis
            at=seq(1, length(actual.freqs$counts)) - 0.5,
            labels=seq(length(actual.freqs$counts), 1, by=-1),
            cex.axis=2
        )
        # Special handling for the y axis of p2.
        if (subject != "p2") {
            axis(2, cex.axis=2)
        } else {
            mark.values=seq(0, 0.3, by=0.1)
            axis(
                2,
                at=mark.values,
                labels=as.character(mark.values),
                cex.axis=2
            )
        }


        # Plot the actual observed counts.
        lines(
            seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
            actual.freqs$counts / total.count,
            type="h",
            lwd=20,
            col="blue",
            lend="butt"
        )

        # Overlay the GLM fitted value and error bars.
        lines(
            rev(x.fit.values$x) - 0.5,
            predictions / total.count,
            lwd=4,
            col="red"
        )
        lines(
            rev(x.fit.values$x) - 0.5,
            upper.bounds / total.count,
            lty=2,
            lwd=2,
            col="red"
        )
        lines(
            rev(x.fit.values$x) - 0.5,
            lower.bounds / total.count,
            lty=2,
            lwd=2,
            col="red"
        )

        # Add text listing the half life.
        decay.rate.summary <- summary(decay.rate.regression)
        x.coef <- decay.rate.summary$coefficients[2, 1]
        x.se <- decay.rate.summary$coefficients[2, 2]
        
        # t_{1/2} = - log(2) / x.coef
        half.life <- - log(2) / x.coef
        half.life.upper <- "\u221e"
        if (x.coef + 1.96 * x.se < 0) {
            half.life.upper <- round(- log(2) / (x.coef + 1.96 * x.se), digits=2)
        }
        half.life.lower <- - log(2) / (x.coef - 1.96 * x.se)

        text(
            x=0,
            y=max.y * 0.90,
            label=substitute(
                paste(
                    t[1/2],
                    " = ",
                    half.life.formatted,
                    " yr",
                    sep=""
                ),
                list(half.life.formatted=round(half.life, digits=2))
            ),
            pos=4,
            cex=2
        )

        text(
            x=0,
            y=max.y * 0.825,
            label=paste0(
                "(95% CI (",
                round(half.life.lower, digits=2),
                ", ",
                half.life.upper,  # we either already rounded it, or it's the infinity symbol
                "))"
            ),
            pos=4,
            cex=2
        )
 
        dev.off()
    }
}

