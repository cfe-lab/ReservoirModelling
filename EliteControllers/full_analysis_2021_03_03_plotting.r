# This is commented out as we will often have already loaded it when we run this script.
# load("full_analysis_2021_03_03.RData")
library(plotrix)

# Make plots of the decay rate estimates.
gaps <- list()
x.tick.marks <- list()
for (subject in subjects) {
    gaps[[subject]] <- c(6500, 49500)
    x.tick.marks[[subject]] <- c(seq(0, 6000, 1000), seq(50000, max(possible.half.lives), 1000))
}

for (subject in subjects) {
    for (regime in regimes) {
        lls <- all.log.likelihoods[[subject]][[regime]]

        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]

        # Filter out values in the plot gap.
        plot.gap <- gaps[[subject]]
        plot.indices <- possible.half.lives < plot.gap[1] | possible.half.lives > plot.gap[2] + 250

        pdf(paste(subject, "_", regime, "_known_vl.pdf", sep=""))
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
            xtics=x.tick.marks[[subject]],
            xticlab=x.tick.marks[[subject]]
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

        if (subject != "p2") {
            # Some special handling if the MLE is bigger than where we put the gap.
            mle.plot <- mle
            if (mle > plot.gap[2]) {
                mle.plot <- mle - (plot.gap[2] - plot.gap[1])
            }
            abline(v=mle.plot, lty="dashed")

            # We get an estimate of the variance of the MLE via the Fisher information.
            estimated.variance <- 
                1 / fisher.information(max.idx, ode.solutions.bin.30[[subject]][[regime]]$bin.freqs)

            lower.bound <- max(0, mle - 2 * sqrt(estimated.variance) * bin.size)
            upper.bound <- mle + 2 * sqrt(estimated.variance) * bin.size
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
                labels="best fit achieved\nwith no decay",
                cex=2
            )
        }

        dev.off()
    }
}

# Plots of reservoir composition.
for (subject in subjects) {
    curr.data <- all.subjects[[subject]]$integration

    days.pre.therapy <- 
        as.numeric(
            all.subjects[[subject]]$art.initiation - all.subjects[[subject]]$infection.date,
            units="days"
        )

    for (regime in regimes) {
        reservoir.dist <- ode.solutions.bin.365[[subject]][[regime]]

        lls <- all.log.likelihoods[[subject]][[regime]]
        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]

        dist.44mo.decay <- decay.distribution(reservoir.dist$bin.freqs, 44 * 30, 365)
        dist.140mo.decay <- decay.distribution(reservoir.dist$bin.freqs, 140 * 30, 365)
        dist.best.fit <- decay.distribution(reservoir.dist$bin.freqs, mle, 365)

        breakpoints = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=365)
        if (!(max(curr.data$days.before.art) %in% breakpoints)) {
            breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + 365)
        }
        actual.freqs <- hist(
            curr.data$days.before.art,
            breaks=breakpoints,
            plot=FALSE
        )
        emp.dist <- actual.freqs$counts / sum(actual.freqs$counts)

        max.y <- max(
            c(
                dist.44mo.decay$bin.dist, 
                dist.140mo.decay$bin.dist, 
                dist.best.fit$bin.dist, 
                reservoir.dist$bin.dist.no.decay,
                emp.dist
            )
        )

        pdf(paste("composition_", subject, "_", regime, "_known_vl.pdf", sep=""))
        # par(mar=c(5, 4, 4, 12) + 0.1, xpd=TRUE)
        par(mar=c(6.5, 6.5, 2, 2) + 0.1)

        # Some defaults, and then some customization for p3.
        y.limits <- c(0, max.y)
        x.label <- "Year prior to ART initiation"
        x.label.cex <- 3
        if (subject == "p3") {
            y.limits[2] <- max.y + 0.2
            x.label <- "Year prior to last viremic episode"
            x.label.cex <- 2.25
        }

        plot(
            c(0, length(emp.dist)),
            y.limits,
            # main=paste(
            #     "Reservoir composition (VL known) for ", 
            #     subject,
            #     " (",
            #     regime,
            #     " case)",
            #     sep=""
            # ),
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

        if (subject != "p2") {
            lines(
                c(0, seq(length(emp.dist) - length(reservoir.dist$bin.dist.no.decay) + 1, length(emp.dist))),
                c(reservoir.dist$bin.dist.no.decay[1], reservoir.dist$bin.dist.no.decay),
                type="S",
                col=colour.alpha.helper("orange", alpha=175, max.colour.value=255),
                lwd=3
            )
        }

        lines(
            c(0, seq(length(emp.dist) - length(dist.140mo.decay$bin.dist) + 1, length(emp.dist))),
            c(dist.140mo.decay$bin.dist[1], dist.140mo.decay$bin.dist),
            type="S",
            col=colour.alpha.helper("grey", alpha=175, max.colour.value=255),
            lwd=3
        )

        lines(
            c(0, seq(length(emp.dist) - length(dist.44mo.decay$bin.dist) + 1, length(emp.dist))),
            c(dist.44mo.decay$bin.dist[1], dist.44mo.decay$bin.dist),
            type="S",
            col="black",
            lty="dotdash",
            lwd=3
        )

        if (subject != "p2") {
            lines(
                c(0, seq(length(emp.dist) - length(dist.best.fit$bin.dist) + 1, length(emp.dist))),
                c(dist.best.fit$bin.dist[1], dist.best.fit$bin.dist),
                type="S",
                col=colour.alpha.helper("green", alpha=175, max.colour.value=255),
                lwd=3
            )
        } else {
            lines(
                c(0, seq(length(emp.dist) - length(reservoir.dist$bin.dist.no.decay) + 1, length(emp.dist))),
                c(reservoir.dist$bin.dist.no.decay[1], reservoir.dist$bin.dist.no.decay),
                type="S",
                col=colour.alpha.helper("green", alpha=175, max.colour.value=255),
                lwd=3
            )
        }

        legend.captions <- c(
            "observed",
            "no decay",
            "140mo decay",
            "44mo decay",
            paste("best-fit decay (", mle, " days)", sep="")
        )
        legend.colours <- c(
            "blue",
            "orange",
            "grey",
            "black",
            "green"
        )
        legend.line.types <- c(
            "solid",
            "solid",
            "solid",
            "dotdash",
            "solid"
        )
        legend.line.widths <- c(20, 3, 3, 3, 3)

        if (subject != "p2") {
            legend(
                "topleft",
                legend=legend.captions,
                col=legend.colours,
                lty=legend.line.types,
                lwd=legend.line.widths,
                cex=1.5,
                bty="n"
            )
        } else {
            # We customize the legend for p2.
            p2.legend.indices <- c(1, 3, 4, 5)
            p2.captions <- legend.captions[p2.legend.indices]
            p2.captions[4] <- "best fit (no decay)"

            legend(
                "topleft",
                legend=p2.captions,
                col=legend.colours[p2.legend.indices],
                lty=legend.line.types[p2.legend.indices],
                lwd=legend.line.widths[p2.legend.indices],
                cex=1.5,
                bty="n"
            )
        }
        dev.off()
    }
}
