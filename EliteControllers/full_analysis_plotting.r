# This is commented out as we will often have already loaded it when we run this script.
# load("full_analysis.RData")
library(plotrix)

# Make plots of the decay rate estimates.
gaps <- list()
x.tick.marks <- list()
for (subject in subjects) {
    gaps[[subject]] <- c(7500, 49500)
    x.tick.marks[[subject]] <- c(seq(0, 7000, 1000), seq(50000, max(possible.half.lives), 1000))
}

bin.size <- 30
mles <- NULL

for (subject in subjects) {
    for (regime in regimes) {
        lls <- all.log.likelihoods[[subject]][[regime]]

        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]
        # We get an estimate of the variance of the MLE via the Fisher information.
        estimated.variance <- 
            1 / fisher.information(max.idx, ode.solutions.bin.30[[subject]][[regime]]$bin.freqs)

        lower.bound <- max(0, mle - 2 * sqrt(estimated.variance) * bin.size)
        upper.bound <- mle + 2 * sqrt(estimated.variance) * bin.size

        mles <- rbind(
            mles,
            data.frame(
                subject=subject,
                regime=regime,
                mle=mle,
                lower.bound=lower.bound,
                upper.bound=upper.bound
            )
        )

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

        mle.row <- mles[mles$subject == subject & mles$regime == regime,]
        mle <- mle.row$mle
        lower.bound <- mle.row$lower.bound
        upper.bound <- mle.row$upper.bound

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

        max.y <- 1

        pdf(paste("composition_", subject, "_", regime, "_known_vl.pdf", sep=""))
        par(mar=c(6.5, 6.5, 2, 2) + 0.1)

        # Some defaults, and then some customization for p3.
        y.limits <- c(0, max.y)
        x.label <- "Year prior to ART initiation"
        x.label.cex <- 2.25
        if (subject == "p3") {
            x.label <- "Year prior to last viremic episode"
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

        comp.line.width <- 6

        if (subject != "p2") {
            lines(
                c(0, seq(length(emp.dist) - length(reservoir.dist$bin.dist.no.decay) + 1, length(emp.dist))),
                c(reservoir.dist$bin.dist.no.decay[1], reservoir.dist$bin.dist.no.decay),
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

        if (subject != "p2") {
            lines(
                c(0, seq(length(emp.dist) - length(dist.best.fit$bin.dist) + 1, length(emp.dist))),
                c(dist.best.fit$bin.dist[1], dist.best.fit$bin.dist),
                type="S",
                col=colour.alpha.helper("green", alpha=175, max.colour.value=255),
                lwd=comp.line.width
            )
        } else {
            lines(
                c(0, seq(length(emp.dist) - length(reservoir.dist$bin.dist.no.decay) + 1, length(emp.dist))),
                c(reservoir.dist$bin.dist.no.decay[1], reservoir.dist$bin.dist.no.decay),
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
        legend.line.widths <- c(20, comp.line.width, comp.line.width, comp.line.width, comp.line.width, 0)
        legend.location <- "topleft"
        if (subject == "p3") {
            legend.location <- "topright"
        }

        if (subject != "p2") {
            legend.coords <- legend(
                legend.location,
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
                legend.location,
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
