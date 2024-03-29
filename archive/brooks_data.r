brooks.data <- read.csv("BrooksData.csv")[,1:8]

names(brooks.data) <- c(
    "pid",
    "est.infection.date",
    "art.start.date",
    "seq.id",
    "collection.date",
    "integration.date.est",
    "integration.date.lower",
    "integration.date.upper"
)

for (col.idx in c(2, 3, 5, 6, 7, 8)) {
    brooks.data[[col.idx]] <- strptime(brooks.data[[col.idx]], "%d-%b-%y")
}

# Get the gap between infection date and ART initiation.
brooks.data$untreated.period <-
    as.numeric(brooks.data$art.start.date - brooks.data$est.infection.date, units="days")

# Get the number of days prior to ART initiation that each sequence is estimated
# to have been introduced to the reservoir.
brooks.data$days.before.art.raw <-
    as.numeric(brooks.data$art.start.date - brooks.data$integration.date.est, units="days")

brooks.data$days.before.art <-
    sapply(
        brooks.data$days.before.art.raw,
        function (x) {
            if (x >= 0) {
                return(x)
            }
            return(0)
        }
    )

integration.data <- brooks.data[, c(1, 5, 9, 11)]

source("reservoir_helpers.r")
library(parallel)

reservoir.dists <- list()
bin.size <- 30
for (i in 1:length(levels(integration.data$pid))) {
    curr.pid <- levels(integration.data$pid)[i]

    curr.pid.data <- integration.data[integration.data$pid == curr.pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    cat("Processing data for PID ", curr.pid, ".\n", sep="")
    cat("Calculating undecayed distribution of reservoir....\n")
    run.time <- system.time(
        reservoir.dist <- undecayed.reservoir.distribution(
            days.pre.therapy,
            bin.size=bin.size
        )
    )
    cat("Time elapsed:\n")
    print(run.time)
    reservoir.dists[[curr.pid]] <- reservoir.dist
}

possible.half.lives <- (1:200) * bin.size  # months, roughly
all.log.likelihoods <- list()  # this will be doubly-indexed by pid and collection date

for (i in 1:length(levels(integration.data$pid))) {
    curr.pid <- levels(integration.data$pid)[i]
    all.log.likelihoods[[curr.pid]] <- list()

    curr.pid.data <- integration.data[integration.data$pid == curr.pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    reservoir.dist <- reservoir.dists[[curr.pid]]

    all.collection.dates <- as.character(unique(curr.pid.data$collection.date))
    for (col.date in all.collection.dates) {
        curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
        cat(
            "Calculating likelihoods on decayed distributions for collection date ",
            col.date,  # remember this has been converted to a string already
            "....\n",
            sep=""
        )
        run.time <- system.time(
            log.likelihoods.no.factorial <- mclapply(
                possible.half.lives,
                function (x) {
                    cat("Calculating distribution of reservoir with a", x, "day half life....\n")
                    decayed <- decay.distribution(reservoir.dist$bin.freqs, x, bin.size)
                    return(log.likelihood.no.factorial(curr.data$days.before.art, decayed$bin.dist, bin.size))
                },
                mc.cores=3
            )
        )
        # Compute the log factorial term for the likelihood.
        all.log.likelihoods[[curr.pid]][[col.date]] <- (
            unlist(log.likelihoods.no.factorial) + 
            sum(log(1:nrow(curr.data)))
        )

        cat(
            "Likelihood calculations for PID ", 
            curr.pid, 
            ", collection date ", 
            col.date, 
            " complete.\n",
            sep=""
        )
        cat("Time elapsed:\n")
        print(run.time)
    }
}

save.image("brooks_data.RData")

for (pid in names(all.log.likelihoods)) {
    for (col.date in names(all.log.likelihoods[[pid]])) {
        lls <- all.log.likelihoods[[pid]][[col.date]]

        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]

        pdf(paste(pid, "_", col.date,  ".pdf", sep=""))
        plot(
            possible.half.lives,
            lls,
            main=paste(
                "Log likelihoods by decay rate: ", pid, ", collection date ", col.date,
                sep=""
            ),
            xlab="Reservoir half life (days)",
            ylab="Log likelihood"
        )
        abline(v=mle, lty="dashed")

        # We get an estimate of the variance of the MLE via the Fisher information.
        estimated.variance <- 1 / fisher.information(max.idx, reservoir.dists[[pid]]$bin.freqs)

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


# Make plots where, for each bin, we show:
# - a barplot showing the number of actual sequences are in that bin
# - a line for a "hypothetical" reservoir with 44mo decay
# - a line for a hypothetical reservoir with 140mo decay
# - a line with the best-fit decay rate
bin.size = 365
reservoir.dists.1yr.bins <- list()
for (i in 1:length(levels(integration.data$pid))) {
    curr.pid <- levels(integration.data$pid)[i]

    curr.pid.data <- integration.data[integration.data$pid == curr.pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    cat("Processing data for PID ", curr.pid, ".\n", sep="")
    cat("Calculating undecayed distribution of reservoir....\n")
    run.time <- system.time(
        reservoir.dist <- undecayed.reservoir.distribution(
            days.pre.therapy,
            bin.size=bin.size
        )
    )
    cat("Time elapsed:\n")
    print(run.time)
    reservoir.dists.1yr.bins[[curr.pid]] <- reservoir.dist
}

for (i in 1:length(levels(integration.data$pid))) {
    curr.pid <- levels(integration.data$pid)[i]
    curr.pid.data <- integration.data[integration.data$pid == curr.pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    reservoir.dist <- reservoir.dists.1yr.bins[[curr.pid]]

    all.collection.dates <- as.character(unique(curr.pid.data$collection.date))
    for (col.date in all.collection.dates) {
        curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]

        lls <- all.log.likelihoods[[curr.pid]][[col.date]]
        max.idx <- which.max(lls)
        mle <- possible.half.lives[max.idx]

        dist.44mo.decay <- decay.distribution(reservoir.dist$bin.freqs, 44 * 30, bin.size)
        dist.140mo.decay <- decay.distribution(reservoir.dist$bin.freqs, 140 * 30, bin.size)
        dist.best.fit <- decay.distribution(reservoir.dist$bin.freqs, mle, bin.size)

        breakpoints = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=bin.size)
        if (!(max(curr.data$days.before.art) %in% breakpoints)) {
            breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + bin.size)
        }
        actual.freqs <- hist(
            curr.data$days.before.art,
            breaks=breakpoints,
            plot=FALSE
        )
        emp.dist <- actual.freqs$counts / sum(actual.freqs$counts)

        max.y <- max(c(dist.44mo.decay$bin.dist, dist.140mo.decay$bin.dist, dist.best.fit$bin.dist, emp.dist))

        pdf(paste("composition_", curr.pid, "_", col.date,  ".pdf", sep=""))
        par(mar=c(5, 4, 4, 12) + 0.1, xpd=TRUE)

        plot(
            c(0, length(emp.dist)),
            c(0, max.y),
            main=paste(
                "Reservoir composition for ", 
                curr.pid,
                "\n(collection date ",
                col.date,
                ")",
                sep=""
            ),
            xlab="Year prior to ART initiation",
            ylab="Proportion",
            type="n",
            xaxt="n"
        )

        axis(
            1,  # this is the x axis
            at=seq(1, length(emp.dist)) - 0.5,
            labels=seq(length(emp.dist), 1, by=-1)
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

        legend(
            "topright",
            legend=c(
                "empirical",
                "140mo decay",
                "44mo decay",
                paste(mle, "day decay", sep="")
            ),
            col=c(
                "blue",
                "grey",
                "black",
                "green"
            ),
            lty=c(
                "solid",
                "solid",
                "dotdash",
                "solid"
            ),
            lwd=c(20, 3, 3, 3),
            inset=c(-0.5, 0)
        )
        dev.off()
    }
}

