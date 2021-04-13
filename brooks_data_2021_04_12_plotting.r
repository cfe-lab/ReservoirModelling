# This builds on the plotting code from 2021-03-25 but with the visuals
# tweaked.

# This is commented out because we will often already have this loaded.
# load("brooks_data_2021_03_25.RData")


# Make plots of the decay rate estimates.
for (vl.info in c("known", "typical")) {
    for (pid in names(all.log.likelihoods[[vl.info]])) {
        for (col.date in names(all.log.likelihoods[[vl.info]][[pid]])) {
            lls <- all.log.likelihoods[[vl.info]][[pid]][[col.date]]

            col.date.str <- "all collection dates"
            if (col.date != "combined") {
                col.date.str <- paste("collected", col.date)
            }

            max.idx <- which.max(lls)
            mle <- possible.half.lives[max.idx]

            pdf(paste(pid, "_", col.date,  "_", vl.info, "_vl.pdf", sep=""))
            plot(
                possible.half.lives,
                lls,
                main=paste(
                    "LLs (using ",
                    vl.info,
                    " VL) by decay rate: ", 
                    pid, 
                    ", ",
                    col.date.str,
                    sep=""
                ),
                xlab="Reservoir half life (days)",
                ylab="Log likelihood"
            )
            abline(v=mle, lty="dashed")

            # We get an estimate of the variance of the MLE via the Fisher information.
            estimated.variance <- (
                1 / fisher.information(
                    max.idx, 
                    ode.solutions.bin.30[[vl.info]][[pid]]$bin.freqs
                )
            )

            lower.bound <- max(0, mle - 2 * sqrt(estimated.variance) * bin.size)
            upper.bound <- mle + 2 * sqrt(estimated.variance) * bin.size
            abline(v=lower.bound, lty="dotted")
            abline(v=upper.bound, lty="dotted")

            text(
                x=mle,
                y=sum(range(lls)) / 2,
                labels=paste(mle, "days"),
                pos=4
            )
            text(
                x=upper.bound,
                y=min(lls) + (max(lls) - min(lls))/ 4,
                labels=paste(round(upper.bound, digits=2), "days"),
                pos=4
            )

            dev.off()
        }
    }
}


# Plots of reservoir composition.
for (vl.info in c("known", "typical")) {
    for (pid in names(all.log.likelihoods[[vl.info]])) {
        curr.pid.data <- integration.data[integration.data$pid == pid,]
        days.pre.therapy <- curr.pid.data$untreated.period[1]

        # Remember that these are strings, and may also be the string "combined",
        # meaning that information from all collection dates was combined.
        all.collection.dates <- names(all.log.likelihoods[[vl.info]][[pid]])
        for (col.date in all.collection.dates) {
            curr.data <- curr.pid.data
            if (col.date != "combined") {
                curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
            }

            reservoir.dist <- ode.solutions.bin.365[[vl.info]][[pid]]

            lls <- all.log.likelihoods[[vl.info]][[pid]][[col.date]]
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

            max.y <- 1

            pdf(paste("composition_", pid, "_", col.date,  "_", vl.info, "_vl.pdf", sep=""))
            par(mar=c(6.5, 6.5, 2, 2) + 0.1)

            plot(
                c(0, length(emp.dist)),
                c(0, max.y),
                # main=paste(
                #     "Reservoir composition (",
                #     vl.info,
                #     " VL) for ", 
                #     curr.pid,
                #     "\n(collected ",
                #     col.date,
                #     ")",
                #     sep=""
                # ),
                xlab=NA,
                ylab=NA,
                type="n",
                xaxt="n",
                cex.axis=2
            )

            title(
                xlab="Year prior to ART initiation",
                line=3.5,
                cex.lab=2.25
            )

            title(
                ylab="Proportion",
                line=3.5,
                cex.lab=3
            )

            axis(
                1,  # this is the x axis
                at=seq(1, length(emp.dist)) - 0.5,
                labels=seq(length(emp.dist), 1, by=-1),
                cex.axis=2
            )

            if (col.date != "combined") {
                lines(
                    seq(length(emp.dist), 1, by=-1) - 0.5,
                    emp.dist,
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
                    seq(length(emp.dist), 1, by=-1) - 0.5,
                    first.freqs$counts / total,
                    type="h",
                    lwd=20,
                    col="steelblue",
                    lend="butt"
                )

                # Then plot the later ones on top.
                segments(
                    x0=seq(length(emp.dist), 1, by=-1) - 0.5,
                    y0=first.freqs$counts / total,
                    y1=(first.freqs$counts + last.freqs$counts) / total,
                    lwd=20,
                    col="turquoise",
                    lend="butt"
                )
            }

            lines(
                c(0, seq(length(emp.dist) - length(dist.140mo.decay$bin.dist) + 1, length(emp.dist))),
                c(dist.140mo.decay$bin.dist[1], dist.140mo.decay$bin.dist),
                type="S",
                col="grey",
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
                col="green",
                lwd=3
            )

            # Build up the legend, with some customization for the "combined" cases.
            legend.labels <- "sampled"
            legend.colours <- "blue"
            legend.ltys <- "solid"
            legend.lwds <- 20
            if (col.date == "combined") {
                legend.labels <- c(
                    "sampling point 1",
                    "sampling point 2"
                )
                legend.colours <- c("steelblue", "turquoise")
                legend.ltys <- c("solid", "solid")
                legend.lwds <- c(20, 20)
            }
            legend.labels <- c(
                legend.labels,
                "140mo decay",
                "44mo decay",
                paste(mle, "day decay", sep="")
            )
            legend.colours <- c(
                legend.colours,
                "grey",
                "black",
                "green"
            )
            legend.ltys <- c(
                legend.ltys,
                "solid",
                "dotdash",
                "solid"
            )
            legend.lwds <- c(legend.lwds, 3, 3, 3)

            legend(
                "topleft",
                legend=legend.labels,
                col=legend.colours,
                lty=legend.ltys,
                lwd=legend.lwds,
                bty="n",
                cex=1.5
            )

            text(
                x=length(emp.dist) * 0.8,
                y=0.95,
                labels=pid,
                cex=2
            )
            dev.off()
        }
    }
}
