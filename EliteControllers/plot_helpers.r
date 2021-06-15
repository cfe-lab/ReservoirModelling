ll.plot <- function(
    lls,
    possible.half.lives,
    mle,
    lower.bound,
    upper.bound,
    plot.gap,
    x.tick.marks
) {
    # Filter out values in the plot gap.
    plot.indices <- possible.half.lives < plot.gap[1] | possible.half.lives > plot.gap[2] + 250

    par(mar=c(6.5, 6.5, 2, 2) + 0.1)
    gap.plot(
        possible.half.lives[plot.indices],
        lls[plot.indices],
        gap=plot.gap,
        gap.axis="x",
        # main=paste("LLs (using known VL) by decay rate: ", subject, " (", regime, " case)", sep=""),
        xlab=NA,
        ylab=NA,
        cex.lab=3,
        cex.axis=2,
        xtics=x.tick.marks,
        xticlab=x.tick.marks
    )
    axis.break(1, plot.gap[1], breakcol="snow", style="gap")
    axis.break(1, plot.gap[1] * (1.02), breakcol="black", style="slash")
    axis.break(3, plot.gap[1] * (1.02), breakcol="black", style="slash")

    title(
        xlab="Reservoir half life (days)",
        ylab="Log likelihood",
        line=3.5,
        cex.lab=3
    )

    # if (subject != "p2" || regime != "min") {
    if (mle != max(possible.half.lives)) {
        # Some special handling if the MLE is bigger than where we put the gap.
        mle.plot <- mle
        if (mle > plot.gap[2]) {
            mle.plot <- mle - (plot.gap[2] - plot.gap[1])
        }
        abline(v=mle.plot, lty="dashed")
        abline(v=lower.bound, lty="dotted")

        # Some special handling if the upper bound is bigger than where we put the gap.
        upper.bound.plot <- upper.bound
        if (upper.bound > plot.gap[2]) {
            upper.bound.plot <- upper.bound - (plot.gap[2] - plot.gap[1])
        }
        abline(v=upper.bound.plot, lty="dotted")

        position <- 4
        text(
            x=mle.plot,
            y=sum(range(lls)) / 2,
            labels=paste(mle, "days"),
            pos=position,
            cex=2,
            offset=0.1
        )

        text(
            x=upper.bound.plot,
            y=min(lls) + (max(lls) - min(lls))/ 4,
            labels=paste(round(upper.bound, digits=2), "days"),
            pos=position,
            cex=2,
            offset=0.1
        )
    } else {
        text(
            x=51000 - (plot.gap[2] - plot.gap[1]),
            y=sum(range(lls)) / 2,
            labels="no meaningful\nMLE found",
            cex=2
        )
    }
}


composition.plot <- function(
    bin.freqs,  # the 365-day-binned bin frequencies
    days.before.art,  # the actual integration dates, in days-before-ART
    days.pre.therapy,  # the number of days the individual spent untreated
    mle,
    lower.bound,
    upper.bound,
    bayes.factor=NA,  # if this isn't NA, plot it in the legend
    max.y=1,
    best.fit.no.decay=FALSE,
    x.label="Year prior to ART initiation",
    x.label.cex=2.25,
    legend.location="topleft",
    comp.line.width=6
) {
    dist.44mo.decay <- decay.distribution(bin.freqs, 44 * 30, 365)
    dist.140mo.decay <- decay.distribution(bin.freqs, 140 * 30, 365)
    dist.best.fit <- decay.distribution(bin.freqs, mle, 365)
    bin.dist.no.decay <- bin.freqs / sum(bin.freqs)

    breakpoints <- compute.bin.breakpoints(
        days.before.art,
        days.pre.therapy,
        bin.size=365
    )
    actual.freqs <- hist(
        days.before.art,
        breaks=breakpoints,
        plot=FALSE
    )
    emp.dist <- actual.freqs$counts / sum(actual.freqs$counts)

    par(mar=c(6.5, 6.5, 2, 2) + 0.1)
    plot(
        c(0, length(emp.dist)),
        c(0, max.y),
        xlab=NA,
        ylab=NA,
        type="n",
        xaxt="n",
        cex.lab=3,
        cex.axis=2
    )

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
        at=seq(1, length(emp.dist)) - 0.5,
        labels=seq(length(emp.dist), 1, by=-1),
        cex.axis=2
    )

    lines(
        seq(length(emp.dist), 1, by=-1) - 0.5,
        emp.dist,
        type="h",
        lwd=20,
        col="blue",
        lend="butt"
    )

    comp.line.width <- 6

    if (!best.fit.no.decay) {
        lines(
            c(0, seq(length(emp.dist) - length(bin.dist.no.decay) + 1, length(emp.dist))),
            c(bin.dist.no.decay[1], bin.dist.no.decay),
            type="S",
            col=colour.alpha.helper("orange", alpha=175, max.colour.value=255),
            lwd=comp.line.width
        )
    }

    lines(
        c(0, seq(length(emp.dist) - length(dist.140mo.decay$bin.dist) + 1, length(emp.dist))),
        c(dist.140mo.decay$bin.dist[1], dist.140mo.decay$bin.dist),
        type="S",
        col=colour.alpha.helper("grey", alpha=175, max.colour.value=255),
        lwd=comp.line.width
    )

    lines(
        c(0, seq(length(emp.dist) - length(dist.44mo.decay$bin.dist) + 1, length(emp.dist))),
        c(dist.44mo.decay$bin.dist[1], dist.44mo.decay$bin.dist),
        type="S",
        col="black",
        lty="dotdash",
        lwd=comp.line.width
    )

    if (!best.fit.no.decay) {
        lines(
            c(0, seq(length(emp.dist) - length(dist.best.fit$bin.dist) + 1, length(emp.dist))),
            c(dist.best.fit$bin.dist[1], dist.best.fit$bin.dist),
            type="S",
            col=colour.alpha.helper("green", alpha=175, max.colour.value=255),
            lwd=comp.line.width
        )
    } else {
        lines(
            c(0, seq(length(emp.dist) - length(bin.dist.no.decay) + 1, length(emp.dist))),
            c(bin.dist.no.decay[1], bin.dist.no.decay),
            type="S",
            col=colour.alpha.helper("green", alpha=175, max.colour.value=255),
            lwd=comp.line.width
        )
    }

    # Convert the MLE from days to years.  Our definition of a "year" is simply
    # 365 days (i.e. we ignored leap years).
    legend.captions <- c(
        "observed",
        "no decay",
        expression(paste(t[1/2], " = 140 mo")),
        expression(paste(t[1/2], " = 44 mo")),
        substitute(
            paste(
                "best-fit ",
                t[1/2],
                " = ",
                best.fit.decay,
                " years",
                sep=""
            ),
            list(
                best.fit.decay=round(mle / 365, digits=2)
            )
        ),
        paste0(
            "(95% CI [", 
            round(lower.bound / 365, digits=2), 
            ", ", 
            round(upper.bound / 365, digits=2), 
            "))"
        )
    )
    legend.colours <- c(
        "blue",
        "orange",
        "grey",
        "black",
        "green",
        "black"  # this should be ignored
    )
    legend.line.types <- c(
        "solid",
        "solid",
        "solid",
        "dotdash",
        "solid",
        NA
    )
    legend.line.widths <- c(20, rep(comp.line.width, 4), 0)

    if (best.fit.no.decay) {
        # We customize the legend for cases where no meaningful MLE was found.
        p2.legend.indices <- c(1, 3, 4, 5)
        legend.captions <- legend.captions[p2.legend.indices]
        legend.captions[4] <- "best fit (no decay)"
        legend.colours <- legend.colours[p2.legend.indices]
        legend.line.types <- legend.line.types[p2.legend.indices]
        legend.line.widths <- legend.line.widths[p2.legend.indices]
    }

    if (!is.na(bayes.factor)) {
        # Add a line to the legend.
        legend.captions <- c(legend.captions, paste0("BF = ", format(bayes.factor, digits=3, scientific=0)))
        legend.colours <- c(legend.colours, "black")  # the last entry should be ignored
        legend.line.types <- c(legend.line.types, NA)
        legend.line.widths <- c(legend.line.widths, 0)
    }

    legend(
        legend.location,
        legend=legend.captions,
        col=legend.colours,
        lty=legend.line.types,
        lwd=legend.line.widths,
        cex=1.5,
        bty="n"
    )
}


model.free.plot <- function(
    decay.rate.regression,
    actual.counts,  # plug in actual.freqs$counts
    max.y=1,
    x.label="Year prior to ART initiation"
) {
    # Make some plots of the regressions.
    # The error bars were computed with the guidance of this helpful blog post: 
    # https://www.r-bloggers.com/2018/12/confidence-intervals-for-glms/
    x.fit.values <- data.frame(x=seq(0.5, length(actual.counts) + 0.5, by=0.1))
    fit <- predict(decay.rate.regression, newdata=x.fit.values, se.fit=TRUE)

    predictions <- exp(fit$fit)
    upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
    lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)

    total.count <- sum(actual.counts)

    par(mar=c(6.5, 6.5, 2, 2) + 0.1)
    plot(
        c(0, length(actual.counts)),
        c(0, max.y),
        xlab=NA,
        ylab=NA,
        type="n",
        xaxt="n",
        yaxt="n",
        cex.lab=3,
        cex.axis=2
    )

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
        at=seq(1, length(actual.counts)) - 0.5,
        labels=seq(length(actual.counts), 1, by=-1),
        cex.axis=2
    )
    axis(2, cex.axis=2)

    # Plot the actual observed counts.
    lines(
        seq(length(actual.counts), 1, by=-1) - 0.5,
        actual.counts / total.count,
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
}
