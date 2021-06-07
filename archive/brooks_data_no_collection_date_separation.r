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

integration.data <- brooks.data[, c(1, 9, 11)]

source("reservoir_helpers.r")
library(parallel)

bin.size <- 30
possible.half.lives <- (1:200) * bin.size  # months, roughly
all.log.likelihoods <- list()

for (i in 1:length(levels(integration.data$pid))) {
    curr.pid <- levels(integration.data$pid)[i]
    curr.data <- integration.data[integration.data$pid == curr.pid,]
    days.pre.therapy <- curr.data$untreated.period[i]

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

    cat("Calculating likelihoods on decayed distributions....\n")
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
    all.log.likelihoods[[curr.pid]] <- unlist(log.likelihoods.no.factorial) + sum(log(1:nrow(curr.data)))

    cat("Likelihood calculations for PID", curr.pid, "complete.\n")
    cat("Time elapsed:\n")
    print(run.time)
}

save.image("brooks_data_no_collection_date_separation.RData")

for (pid in names(all.log.likelihoods)) {
    lls <- all.log.likelihoods[[pid]]

    max.idx <- which.max(lls)
    mle <- possible.half.lives[max.idx]

    pdf(paste(pid, ".pdf", sep=""))
    plot(
        possible.half.lives,
        lls,
        main=paste("Log likelihoods for subject", pid, "by decay rate"),
        xlab="Reservoir half life (days)",
        ylab="Log likelihood"
    )
    abline(v=mle, lty="dashed")

    text(
        x=mle,
        y=sum(range(lls)) / 2,
        labels=paste(mle, "days"),
        pos=4
    )

    dev.off()
}
