source("../reservoir_helpers.r")


# This helper loads the viral load data.
prepare.vl.data <- function(
    subjects,
    p3.special.handling=TRUE
) {
    all.subjects <- list()

    # Some special dates that we will have to remember for p3.
    p3.art.initiation <- NULL
    p3.blip <- NULL

    for (subject in subjects) {
        # The VL CSV file has columns (date, viral.load, comments).
        vl.csv <- read.csv(paste("../../data/EliteControllers_2020_11_30/", subject, "/VL-VL.csv", sep=""))
        vl.csv$date <- strptime(vl.csv$date, "%d/%m/%Y")

        peak.csv <- read.csv(paste("../../data/EliteControllers_2020_11_30/", subject, "/Peak-Peak.csv", sep=""))
        infection.date <- strptime(peak.csv$est.infection.date[1], "%d/%m/%Y")
        art.initiation <- strptime(peak.csv$art.initiation[1], "%d/%m/%Y")

        if (p3.special.handling && subject == "p3") {
            # Implement special handling for p3: we will "start the model" from the
            # time of their last viremic episode (i.e. blip) rather than from their
            # cART initiation.  We record those dates separately for later use.
            p3.art.initiation <- art.initiation
            viremic <- which(vl.csv$comments != " Undetectable")
            p3.blip <- vl.csv$date[viremic[length(viremic)]]
            art.initiation <- p3.blip

            undetectable.before.blip <- (
                vl.csv$comments == " Undetectable" &
                vl.csv$date < p3.blip
            ) 
            vl.csv$viral.load[undetectable.before.blip] <- 0
        }

        all.subjects[[subject]] <- list(
            vl=vl.csv,
            infection.date=infection.date,
            art.initiation=art.initiation
        )
    }

    return(
        list(
            all.subjects=all.subjects,
            p3.art.initiation=p3.art.initiation,
            p3.blip=p3.blip
        )
    )
}


compute.days.before.art <- function(
    integration.date,
    model.start.date,  # usually ART initiation, but could be a viremic episode
    boundaries=NULL  # shunt values between these boundaries to the closest boundary
) {
    days.before.art.raw <-
        as.numeric(model.start.date - integration.date, units="days")
    
    days.before.art <- sapply(
        days.before.art.raw,
        function (x) {
            if (x >= 0) {
                return(x)
            }
            return(0)
        }
    )

    if (!is.null(boundaries) && length(boundaries) >= 2 && boundaries[1] < boundaries[2]) {
        left.boundary <- boundaries[1]
        right.boundary <- boundaries[2]
        days.before.art <- sapply(
            1:length(days.before.art),
            function (idx) {
                int.date <- integration.date[idx]
                if (int.date <= left.boundary || int.date >= right.boundary) {
                    return(days.before.art[idx])
                }

                days.after.left <- as.numeric(int.date - left.boundary)
                days.before.right <- as.numeric(right.boundary - int.date)

                if (days.after.left <= days.before.right) {
                    return(as.numeric(right.bondary - left.boundary))
                } else {
                    return(0)
                }
            }
        )
    }
}


graft.acute.phase <- function(
    all.subjects,
    acute.phase,
    acute.setpoint=1000
) {
    subjects <- names(all.subjects)
    regimes <- names(acute.phase)
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
            acute.vls <- c(0, peak.vl, acute.setpoint)
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
    return(grafted.vls)
}


solve.odes <- function(all.subjects, grafted.vls) {
    require(deSolve)
    subjects <- names(grafted.vls)
    for (subject in subjects) {
        curr <- all.subjects[[subject]]
        infection.date <- curr[["infection.date"]]
        art.initiation <- curr[["art.initiation"]]

        ode.solutions.bin.30[[subject]] <- list()
        ode.solutions.bin.365[[subject]] <- list()

        regimes <- names(grafted.vls[[subject]])
        # Compute the ODEs for this individual *in each regime*, and bin the reservoir into 30 day intervals.
        for (regime in regimes) {
            cat("Solving ODEs for ", subject, " (", regime, " case).\n", sep="")
            ode.solutions.bin.30[[subject]][[regime]] <- undecayed.reservoir.distribution.given.vl(
                grafted.vls[[subject]][[regime]]$time,
                grafted.vls[[subject]][[regime]]$viral.load / 1000,  # convert c/mL to c/uL
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

    return(
        list(
            bin.30=ode.solutions.bin.30,
            bin.365=ode.solutions.bin.365
        )
    )
}


compute.lls <- function(
    all.subjects,
    ode.solutions.bin.30,
    bin.size=30,  # months, roughly
    possible.half.lives=(1:1800) * 30,
    mc.cores=3
) {
    # Find the decay rate that maximizes the likelihood for each individual.
    require(parallel)
    all.log.likelihoods <- list()  # this will be triply-indexed by pid, regime, and collection date

    subjects <- names(all.subjects)
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
                    mc.cores=mc.cores
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

    return(all.log.likelihoods)
}
