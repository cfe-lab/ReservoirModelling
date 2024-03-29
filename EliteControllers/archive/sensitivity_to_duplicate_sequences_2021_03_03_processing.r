# March 3, 2021: testing sensitivity to duplicate proviral sequences.  These sequences
# have their integration dates truncated to the time of cART initiation; in particular
# p3 does *not* have these dates truncated to the time of their last blip.

library(deSolve)

# These constants are taken from the literature.
acute.phase <- list()
acute.phase[["Miura"]] <- list(
    days.to.peak=31,
    peak.vl=53300,
    days.to.undetectable=41
)
regimes <- "Miura"

# Read in the data.
all.subjects <- list()
subjects <- c("p1", "p2", "p3", "p4")

for (subject in subjects) {
    # The VL CSV file has columns (date, viral.load, comments).
    vl.csv <- read.csv(paste("../../data/EliteControllers_2020_11_30/", subject, "/VL-VL.csv", sep=""))
    vl.csv$date <- strptime(vl.csv$date, "%d/%m/%Y")

    peak.csv <- read.csv(paste("../../data/EliteControllers_2020_11_30/", subject, "/Peak-Peak.csv", sep=""))
    infection.date <- strptime(peak.csv$est.infection.date[1], "%d/%m/%Y")
    art.initiation <- strptime(peak.csv$art.initiation[1], "%d/%m/%Y")

    all.subjects[[subject]] <- list(
        vl=vl.csv,
        infection.date=infection.date,
        art.initiation=art.initiation
    )
}


grafted.vls <- list()
for (subject in subjects) {
    vl <- all.subjects[[subject]]$vl
    infection.date <- all.subjects[[subject]]$infection.date

    # Convert all dates to days relative to the (estimated) infection date.
    time <- as.numeric(vl$date - infection.date, units="days")

    # Check: if the last one of these entries is actually *after* our first data point,
    # then discard it.  (We can get away without checking the other dates, as they all
    # happen to be fine for all four subjects.)
    vl.converted.time <- data.frame(
        time=time,
        viral.load=vl$viral.load,
        comments=as.character(vl$comments)
    )

    grafted.vls[[subject]] <- list()
    for (regime in regimes) {
        days.to.peak <- acute.phase[[regime]]$days.to.peak
        days.to.undetectable <- acute.phase[[regime]]$days.to.undetectable
        peak.vl <- acute.phase[[regime]]$peak.vl

        acute.time <- c(
            0,
            days.to.peak,
            days.to.peak + days.to.undetectable
        )
        acute.vls <- c(0, peak.vl, 1000)
        acute.comments <- rep("graft", 3)
        graft.df <- data.frame(
            time=acute.time,
            viral.load=acute.vls,
            comments=acute.comments
        )

        grafted.vl <- NULL
        if (acute.time[1] > time[1]) {
            grafted.vl <- vl.converted.time
        } else if (acute.time[2] > time[1]) {
            grafted.vl <- rbind(
                graft.df[1,],
                vl.converted.time
            )
        } else if (acute.time[3] > time[1]) {
            grafted.vl <- rbind(
                graft.df[1:2,],
                vl.converted.time
            )
        } else {
            grafted.vl <- rbind(
                graft.df,
                vl.converted.time
            )
        }
        grafted.vls[[subject]][[regime]] <- grafted.vl
    }
}


# Read in the sample time data.
for (subject in subjects) {
    # This is where the proviral integration data *with duplicates* lives.
    sample.times.base.path <- "../../data/EliteControllers_2020_11_30"

    # The sample times CSV file has 4 useful columns, toss the rest.
    subject.data <- read.csv(paste(sample.times.base.path, subject, "sample_times.csv", sep="/"))
    subject.data <- subject.data[, c(1, 4, 5, 6)]
    names(subject.data) <- c(
        "id",
        "integration.date.est",
        "integration.date.lower",
        "integration.date.upper"
    )

    for (col.idx in 2:4) {
        subject.data[[col.idx]] <- strptime(subject.data[[col.idx]], "%Y-%m-%d")
    }

    days.before.art.raw <-
        as.numeric(all.subjects[[subject]]$art.initiation - subject.data$integration.date.est, units="days")
    
    subject.data$days.before.art <-
        sapply(
            days.before.art.raw,
            function (x) {
                if (x >= 0) {
                    return(x)
                }
                return(0)
            }
        )

    all.subjects[[subject]][["integration"]] <- subject.data
}


# Solve the ODEs numerically (this is the slow part).
source("../reservoir_helpers.r")
ode.solutions.bin.30 <- list()
ode.solutions.bin.365 <- list()
for (subject in subjects) {
    curr <- all.subjects[[subject]]
    infection.date <- curr[["infection.date"]]
    art.initiation <- curr[["art.initiation"]]

    ode.solutions.bin.30[[subject]] <- list()
    ode.solutions.bin.365[[subject]] <- list()
    # Compute the ODEs for this individual *in each regime*, and bin the reservoir into 30 day intervals.
    for (regime in regimes) {
        cat("Solving ODEs for ", subject, " (", regime, " case).\n", sep="")
        ode.solutions.bin.30[[subject]][[regime]] <- undecayed.reservoir.distribution.given.vl(
            grafted.vls[[subject]][[regime]]$time,
            grafted.vls[[subject]][[regime]]$viral.load,
            infection.date,
            art.initiation,
            bin.size=30,
            grid.size=0.001
        )

        # Now get the same information binned by 365 days.
        ode.df <- ode.solutions.bin.30[[subject]][[regime]]$solution
        latent.by.bin.365 <- bin.helper(ode.df, bin.size=365, grid.size=0.001)
        ode.solutions.bin.365[[subject]][[regime]] <- list(
            solution=ode.df,
            bin.freqs=latent.by.bin.365,
            bin.dist.no.decay=latent.by.bin.365 / sum(latent.by.bin.365)
        )
    }
}


# Find the decay rate that maximizes the likelihood for each individual.
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


save.image("sensitivity_to_duplicate_sequences_2021_03_03.RData")
