# Jan 21, 2021: this is the same analysis as the 2020_12_22 script
# but the plots are tweaked for inclusion in a manuscript.

load("elite_controllers_known_vl_2020_12_22.RData")

# Find the decay rate that maximizes the likelihood for each individual.
# NEW FOR THIS ANALYSIS:
# look up to 1800 months, which is roughly 150 years.
library(parallel)
bin.size <- 30
possible.half.lives <- (1:1800) * bin.size  # months, roughly
all.log.likelihoods <- list()  # this will be triply-indexed by pid, regime, and collection date

for (subject in subjects) {
    all.log.likelihoods[[subject]] <- list()
    curr.data <- all.subjects[[subject]]$integration

    for (regime in regimes) {
        all.log.likelihoods[[subject]][[regime]] <- list()

        reservoir.dist <- ode.solutions.bin.30[[subject]][[regime]]
        cat(
            "Calculating likelihoods on decayed distributions for subject ",
            subject,
            " (",
            regime,
            " case)....\n",
            sep=""
        )
        run.time <- system.time(
            log.likelihoods.no.factorial <- mclapply(
                possible.half.lives,
                function (x) {
                    decayed <- decay.distribution(reservoir.dist$bin.freqs, x, bin.size)
                    return(
                        log.likelihood.no.factorial(
                            curr.data$days.before.art, 
                            decayed$bin.dist, 
                            bin.size
                        )
                    )
                },
                mc.cores=3
            )
        )
        # Compute the log factorial term for the likelihood.
        all.log.likelihoods[[subject]][[regime]] <- (
            unlist(log.likelihoods.no.factorial) + 
            sum(log(1:nrow(curr.data)))
        )

        cat(
            "Likelihood calculations for subject ", 
            subject, 
            " (",
            regime,
            " case) complete.\n",
            sep=""
        )
        cat("Time elapsed:\n")
        print(run.time)
    }
}


# Make plots of the decay rate estimates.
for (subject in subjects) {
    for (regime in regimes) {
        lls <- all.log.likelihoods[[subject]][[regime]]

        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]

        pdf(paste(subject, "_", regime, "_known_vl.pdf", sep=""))
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


# # For p2, we compute a likelihood ratio/Bayes factor to compare the models with
# # - 6000 day decay
# # - no decay at all
# subject <- "p2"
# curr.data <- all.subjects[[subject]]$integration
# reservoir.dist <- ode.solutions.bin.30[[subject]]

# no.decay.ll.no.factorial <- log.likelihood.no.factorial(
#     curr.data$days.before.art,
#     reservoir.dist$bin.dist.no.decay,
#     30
# )
# no.decay.ll <- no.decay.ll.no.factorial + sum(log(1:nrow(curr.data)))

# # The maximizer, which was 6000 days.
# lls <- all.log.likelihoods[[subject]]
# max.idx <- which.max(lls)
# max.ll <- lls[max.idx]

# bayes.factor <- exp(no.decay.ll - max.ll)


# # Read in data from the idealized VL analysis.
# load("idealized_results.RData")
# curr.data.idealized <- all.subjects.idealized[[subject]]$integration
# reservoir.dist.idealized <- ode.solutions.bin.30.idealized[[subject]]

# lls.idealized <- all.log.likelihoods.idealized[[subject]]
# max.idx.idealized <- which.max(lls.idealized)
# max.ll.idealized <- lls.idealized[max.idx.idealized]

# bf.idealized.vs.known <- exp(max.ll - max.ll.idealized)

# no.decay.ll.idealized <- log.likelihood.no.factorial(
#     curr.data.idealized$days.before.art,
#     reservoir.dist.idealized$bin.dist.no.decay,
#     30
# ) + sum(log(1:nrow(curr.data)))
