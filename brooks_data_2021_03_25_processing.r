brooks.data <- read.csv("../data/BrooksData.csv")[,1:8]

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


# Next, read in the VL data.
vl.data <- read.csv("../data/BrooksDataVLs.csv")
vl.data <- vl.data[, 1:6]
names(vl.data) <- c(
    "pid",
    "est.infection.date",
    "vl.date",
    "vl",
    "art.start.date",
    "is.estimated.peak"
)

for (col.idx in c(2, 3, 5)) {
    vl.data[[col.idx]] <- strptime(vl.data[[col.idx]], "%d-%b-%y")
}

# For the VL data, we want to convert all of our dates to days after infection.
# Note that we don't have a *zero* day yet!
vl.data$days.after.infection <- as.numeric(vl.data$vl.date - vl.data$est.infection.date, units="days")


# Solve the ODEs numerically (this is the slow part).
source("reservoir_helpers.r")
ode.known.bin.30 <- list()
ode.known.bin.365 <- list()
ode.typical.bin.30 <- list()
ode.typical.bin.365 <- list()
for (pid in unique(brooks.data$pid)) {
    pid.sample.data <- brooks.data[brooks.data$pid == pid,]
    pid.vl.data <- vl.data[vl.data$pid == pid,]

    infection.date <- pid.vl.data$est.infection.date[1]
    art.initiation <- pid.vl.data$art.start.date[1]

    # Sanity check: for each PID, check that the estimated infection date and cART initiation date match
    # in both data frames.
    est.infection.date.matches <- pid.sample.data$est.infection.date[1] == infection.date
    if (!est.infection.date.matches) {
        cat("PID", pid, "shows inconsistent infection dates!\n")
    }
    art.initiation.matches <- pid.sample.data$art.start.date[1] == art.initiation
    if (!art.initiation.matches) {
        cat("PID", pid, "shows inconsistent cART initiation dates!")
    }

    # Compute the ODE for this individual, and bin the reservoir into 30 day intervals.
    cat("Solving ODEs for ", pid, ".\n", sep="")
    ode.known.bin.30[[pid]] <- undecayed.reservoir.distribution.given.vl(
        c(0, pid.vl.data$days.after.infection),
        c(0, pid.vl.data$vl),
        infection.date,
        art.initiation,
        bin.size=30,
        grid.size=0.001
    )

    # Now get the same information binned by 365 days.
    ode.df <- ode.known.bin.30[[pid]]$solution
    latent.by.bin.365 <- bin.helper(ode.df, bin.size=365, grid.size=0.001)
    ode.known.bin.365[[pid]] <- list(
        solution=ode.df,
        bin.freqs=latent.by.bin.365,
        bin.dist.no.decay=latent.by.bin.365 / sum(latent.by.bin.365)
    )

    # Now do the same with the "typical" VLs.
    cat("Calculating undecayed distribution of reservoir....\n")
    run.time <- system.time(
        ode.typical.bin.30[[pid]] <- undecayed.reservoir.distribution(
            as.numeric(art.initiation - infection.date, units="days"),
            bin.size=30
        )
    )
    cat("Time elapsed:\n")
    print(run.time)

    # Bin this by 365 days also.
    ode.typical.df <- ode.typical.bin.30[[pid]]$solution
    latent.by.bin.typical.365 <- bin.helper(ode.typical.df, bin.size=365, grid.size=0.001)
    ode.typical.bin.365[[pid]] <- list(
        solution=ode.typical.df,
        bin.freqs=latent.by.bin.typical.365,
        bin.dist.no.decay=latent.by.bin.typical.365 / sum(latent.by.bin.typical.365)
    )
}
ode.solutions.bin.30 <- list(
    known=ode.known.bin.30,
    typical=ode.typical.bin.30
)
ode.solutions.bin.365 <- list(
    known=ode.known.bin.365,
    typical=ode.typical.bin.365
)
# save.image("brooks_data_2021_03_25_partial.RData")


# Next, estimate the decay rate via MLE.  We handle the individuals with 
# more than one collection date specially.

# Find the decay rate that maximizes the likelihood for each individual.
library(parallel)
bin.size <- 30
possible.half.lives <- (1:200) * bin.size  # months, roughly

# This will be triply-indexed by:
# VL info: "known" or "typical"
# PID
# collection date (this could be "combined" also)
all.log.likelihoods <- list()  # this will be triply-indexed by pid and collection date

for (vl.info in c("known", "typical")) {
    all.log.likelihoods[[vl.info]] <- list()

    for (curr.pid in unique(integration.data$pid)) {
        all.log.likelihoods[[vl.info]][[curr.pid]] <- list()

        curr.pid.data <- integration.data[integration.data$pid == curr.pid,]

        reservoir.dist <- ode.solutions.bin.30[[vl.info]][[curr.pid]]

        all.collection.dates <- as.character(unique(curr.pid.data$collection.date))
        if (length(all.collection.dates) > 1) {
            all.collection.dates <- c(all.collection.dates, "combined")
        }
        for (col.date in all.collection.dates) {
            curr.data <- curr.pid.data
            col.date.str <- "with collection dates combined"
            if (col.date != "combined") {
                curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
                col.date.str <- paste("for collection date", col.date)
            }
            cat(
                "Calculating likelihoods on decayed distributions ",
                col.date.str,
                "....\n",
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
            all.log.likelihoods[[vl.info]][[curr.pid]][[col.date]] <- (
                unlist(log.likelihoods.no.factorial) + 
                sum(log(1:nrow(curr.data)))
            )

            cat(
                "Likelihood calculations for PID ", 
                curr.pid, 
                " ",
                col.date.str,
                " (",
                vl.info,
                " VL data) complete.\n",
                sep=""
            )
            cat("Time elapsed:\n")
            print(run.time)
        }
    }
}

save.image("brooks_data_2021_03_25.RData")
