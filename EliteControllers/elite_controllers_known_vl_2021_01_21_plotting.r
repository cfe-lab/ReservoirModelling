load("elite_controllers_known_vl_2021_01_21.RData")


# Make plots of the decay rate estimates.
for (subject in subjects) {
    for (regime in regimes) {
        lls <- all.log.likelihoods[[subject]][[regime]]

        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]

        pdf(paste(subject, "_", regime, "_known_vl.pdf", sep=""))
        par(mar=c(8,7,2,2) + 0.1)
        plot(
            possible.half.lives,
            lls,
            # main=paste("LLs (using known VL) by decay rate: ", subject, " (", regime, " case)", sep=""),
            xlab="Reservoir half life (days)",
            ylab="Log likelihood",
            cex.lab=3,
            cex.axis=2
        )
        abline(v=mle, lty="dashed")

        # We get an estimate of the variance of the MLE via the Fisher information.
        estimated.variance <- 
            1 / fisher.information(max.idx, ode.solutions.bin.30[[subject]][[regime]]$bin.freqs)

        lower.bound <- max(0, mle - 2 * sqrt(estimated.variance) * bin.size)
        upper.bound <- mle + 2 * sqrt(estimated.variance) * bin.size
        abline(v=lower.bound, lty="dotted")
        abline(v=upper.bound, lty="dotted")

        text(
            x=mle,
            y=sum(range(lls)) / 2,
            labels=paste(mle, "days"),
            pos=4,
            cex=2
        )
        text(
            x=upper.bound,
            y=min(lls) + (max(lls) - min(lls))/ 4,
            labels=paste(round(upper.bound, digits=2), "days"),
            pos=4,
            cex=2
        )

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
        par(mar=c(8,7,2,2) + 0.1)

        plot(
            c(0, length(emp.dist)),
            c(0, max.y),
            # main=paste(
            #     "Reservoir composition (VL known) for ", 
            #     subject,
            #     " (",
            #     regime,
            #     " case)",
            #     sep=""
            # ),
            xlab="Year prior to ART initiation",
            ylab="Proportion",
            type="n",
            xaxt="n",
            cex.lab=3,
            cex.axis=2
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
            lend="square"
        )

        lines(
            c(0, seq(length(emp.dist) - length(reservoir.dist$bin.dist.no.decay) + 1, length(emp.dist))),
            c(reservoir.dist$bin.dist.no.decay[1], reservoir.dist$bin.dist.no.decay),
            type="S",
            col="orange",
            lwd=3
        )

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

        lines(
            c(0, seq(length(emp.dist) - length(dist.best.fit$bin.dist) + 1, length(emp.dist))),
            c(dist.best.fit$bin.dist[1], dist.best.fit$bin.dist),
            type="S",
            col=colour.alpha.helper("green", alpha=175, max.colour.value=255),
            lwd=3
        )

        legend(
            "topleft",
            legend=c(
                "observed",
                "no decay",
                "140mo decay",
                "44mo decay",
                paste("best-fit decay rate (", mle, ")", sep="")
            ),
            col=c(
                "blue",
                "orange",
                "grey",
                "black",
                "green"
            ),
            lty=c(
                "solid",
                "solid",
                "solid",
                "dotdash",
                "solid"
            ),
            lwd=c(20, 3, 3, 3, 3)
        )
        dev.off()
    }
}
