# This is a refinement of the 2021-03-25 script, with visual tweaks.

# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.

# This is commented out as we will often have already loaded it when we run this script.
# load("brooks_data_2021_03_25.RData")


# This will be doubly-indexed by:
#  - PID
#  - collection date (could also be "combined")
# As the VL data doesn't enter into this analysis, there's nothing to separate the
# "known" and "typical" VL cases.
all.regressions <- list()

vl.info <- "known"
all.regressions <- list()
for (pid in names(all.log.likelihoods[[vl.info]])) {
    all.regressions[[pid]] <- list()

    curr.pid.data <- integration.data[integration.data$pid == pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    for (col.date in names(all.log.likelihoods[[vl.info]][[pid]])) {
        curr.data <- curr.pid.data
        if (col.date != "combined") {
            curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
        }

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
        all.regressions[[vl.info]][[pid]][[col.date]] <- decay.rate.regression

        x.fit.values <- data.frame(x=seq(0.5, length(actual.freqs$counts) + 0.5, by=0.1))
        fit <- predict(decay.rate.regression, newdata=x.fit.values, se.fit=TRUE)

        predictions <- exp(fit$fit)
        upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
        lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)

        total.count <- sum(actual.freqs$counts)
        max.y <- max(c(actual.freqs$counts, upper.bounds)) / total.count
        max.y <- min(1, max.y)

        cairo_pdf(
            paste0(
                "model_free_decay_rate_", 
                pid,
                "_",
                col.date,
                ".pdf"
            )
        )
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
        axis(2, cex.axis=2)

        if (col.date != "combined") {
            lines(
                seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
                actual.freqs$counts / total.count,
                type="h",
                lwd=20,
                col="blue",
                lend="butt"
            )
        } else {
            # Plot the empirical distribution as a stacked bar graph.  
            # There are at most two collection dates.
            collection.dates <- unique(curr.data$collection.date)
            first.collected <- curr.data[curr.data$collection.date == min(collection.dates),]
            last.collected <- curr.data[curr.data$collection.date == max(collection.dates),]

            first.freqs <- hist(
                first.collected$days.before.art,
                breaks=breakpoints,
                plot=FALSE
            )
            last.freqs <- hist(
                last.collected$days.before.art,
                breaks=breakpoints,
                plot=FALSE
            )

            # Sanity check: the sum of these frequencies should equal actual.freqs$count.
            if (any(actual.freqs$counts != first.freqs$counts + last.freqs$counts)) {
                cat(
                    "Warning: counts don't add up for case vl.info=",
                    vl.info,
                    ", pid=",
                    pid,
                    "\n",
                    sep=""
                )
            }

            total <- sum(first.freqs$counts, last.freqs$counts)

            # First plot the earlier ones as with the other plots.
            lines(
                seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
                first.freqs$counts / total,
                type="h",
                lwd=20,
                col="steelblue",
                lend="butt"
            )

            # Then plot the later ones on top.
            segments(
                seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
                y0=first.freqs$counts / total,
                y1=(first.freqs$counts + last.freqs$counts) / total,
                lwd=20,
                col="turquoise",
                lend="butt"
            )
        }

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
            y=max.y * 0.85,
            label=paste(
                pid,
                "\n",
                "Half-life ",
                round(half.life, digits=2),
                " years\n(95% CI (",
                round(half.life.lower, digits=2),
                ", ",
                half.life.upper,  # we either already rounded it, or it's the infinity symbol
                "))",
                sep=""
            ),
            pos=4,
            cex=2
        )

        dev.off()
    }
}
