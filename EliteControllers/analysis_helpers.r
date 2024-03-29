source("../reservoir_helpers.r")


# These constants are taken from the literature.
acute.phase <- list()
acute.phase[["min"]] <- list(
    days.to.peak=14,
    peak.vl=1954,
    days.to.undetectable=6 * 7
)
acute.phase[["median"]] <- list(
    days.to.peak=23,
    peak.vl=12902,
    days.to.undetectable=22 * 7
)
acute.phase[["max"]] <- list(
    days.to.peak=30,
    peak.vl=71550,
    days.to.undetectable=43 * 7
)
acute.phase[["Miura"]] <- list(
    days.to.peak=31,
    peak.vl=53300,
    days.to.undetectable=41
)
regimes <- names(acute.phase)
subjects <- c("p1", "p2", "p3", "p4")


# This helper loads the viral load data.
prepare.vl.data <- function(
    subjects,
    p3.special.handling=TRUE
) {
    vl.data <- list()

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

        vl.data[[subject]] <- list(
            vl=vl.csv,
            infection.date=infection.date,
            art.initiation=art.initiation
        )
    }

    return(
        list(
            vl.data=vl.data,
            p3.art.initiation=p3.art.initiation,
            p3.blip=p3.blip
        )
    )
}


prepare.integration.data.helper <- function(
    subjects,  # which subjects to prepare data for
    paths,  # a list keyed by subject with the paths to the files to use
    vl.data,
    date.columns,  # indices *before reducing to useful columns* of columns which should be converted to dates
    useful.column.indices,  # which columns to retain
    useful.column.names,  # what to rename these remaining columns
    p3.boundaries=NULL  # specify this if you want special handling of P3's VL data
) {
    integration.data <- list()
    for (subject in subjects) {
        subject.data <- read.csv(paths[[subject]])

        for (col.idx in date.columns) {
            subject.data[[col.idx]] <- strptime(subject.data[[col.idx]], "%Y-%m-%d")
        }

        subject.data <- subject.data[, useful.column.indices]
        names(subject.data) <- useful.column.names

        if (subject != "p3" || is.null(p3.boundaries)) {
            subject.data$days.before.art <- compute.days.before.art(
                subject.data$integration.date.est,
                vl.data[[subject]]$art.initiation,
                vl.data[[subject]]$infection.date
            )
        } else {
            subject.data$days.before.art <- compute.days.before.art(
                subject.data$integration.date.est,
                vl.data[[subject]]$art.initiation,
                vl.data[[subject]]$infection.date,
                boundaries=p3.boundaries
            )
        }
        integration.data[[subject]] <- subject.data
    }
    return(integration.data)
}


prepare.integration.data <- function(
    subjects,
    vl.data,
    remove.duplicates,
    p3.boundaries
) {
    paths <- list()
    for (subject in subjects) {
        paths[[subject]] <- paste(
            "../../data/IntegrationData_2021_06_10",
            paste0(
                subject,
                "_integration.csv"
            ),
            sep="/"
        )
    }

    integration.data <- prepare.integration.data.helper(
        subjects,
        paths,
        vl.data,
        4:6,
        c(1, 4, 5, 6, 7),
        c(
            "id",
            "integration.date.est",
            "integration.date.lower",
            "integration.date.upper",
            "duplicate"
        ),
        p3.boundaries
    )

    if (remove.duplicates) {
        for (subject in subjects) {
            subject.data <- integration.data[[subject]]
            integration.data[[subject]] <- subject.data[is.na(subject.data$duplicate),]
        }
    }

    return(integration.data)
}


prepare.alternative.trees.integration.data <- function(subjects, vl.data, p3.boundaries) {
    paths <- list()
    for (subject in subjects) {
        paths[[subject]] <- paste(
            "../../data/Alternative trees_proviral integration dates_21Jul2021",
            paste0(
                subject,
                "_alt_reservoir integration dates_21Jul2021.csv"
            ),
            sep="/"
        )
    }

    return(
        prepare.integration.data.helper(
            subjects,
            paths,
            vl.data,
            4:6,
            c(1, 4, 5, 6),
            c(
                "id",
                "integration.date.est",
                "integration.date.lower",
                "integration.date.upper"
            ),
            p3.boundaries
        )
    )
}


prepare.lsd.integration.data <- function(subjects, vl.data, p3.boundaries) {
    paths <- list()
    data.dir <- "../../data/LSD_Proviral integration dates_21Jul2021"
    for (subject in subjects) {
        if (subject %in% c("p1", "p2", "p3")) {
            paths[[subject]] <- paste(
                data.dir,
                paste0(
                    subject,
                    "_LSD_proviral integration dates_21Jul2021.csv"
                ),
                sep="/"
            )
        } else {  # subject == "p4"
            paths[[subject]] <- paste(
                data.dir,
                paste0(
                    subject,
                    "_LSD_proviral integration dates_22Jul2021.csv"
                ),
                sep="/"
            )
        }
    }

    return(
        prepare.integration.data.helper(
            subjects,
            paths,
            vl.data,
            4:6,
            c(1, 4, 5, 6),
            c(
                "id",
                "integration.date.est",
                "integration.date.lower",
                "integration.date.upper"
            ),
            p3.boundaries
        )
    )
}


prepare.nearest.neighbour.integration.data <- function(subjects, vl.data) {
    paths <- list()
    for (subject in subjects) {
        paths[[subject]] <- paste(
            "../../data/NearestNeighbourEstimation_2021_06_08", 
            paste0(
                subject,
                "_integration.csv"
            ),
            sep="/"
        )
    }

    # These files have two useful columns (the first and second):
    # PBMC ID
    # NN Year (which is actually a date in YYYY-MM-DD format)
    return(
        prepare.integration.data.helper(
            subjects,
            paths,
            vl.data,
            2,  # convert the 2nd column to a date
            1:2,  # retain only the first two columns
            c("id", "integration.date.est"),
            NULL  # no special handling of P3's integration dates
        )
    )
}


compute.days.before.art <- function(
    integration.date,
    model.start.date,  # usually ART initiation, but could be a viremic episode
    earliest.possible.date=NULL,  # usually EDI
    boundaries=NULL  # shunt values between these boundaries to the closest boundary
) {
    days.before.art.raw <-
        as.numeric(model.start.date - integration.date, units="days")

    earliest.possible.days.before.art <- Inf
    if (!is.null(earliest.possible.date)) {
        earliest.possible.days.before.art <-
            as.numeric(model.start.date - earliest.possible.date, units="days")
    }
    
    days.before.art <- sapply(
        days.before.art.raw,
        function (x) {
            if (x >= 0) {
                if (x > earliest.possible.days.before.art) {
                    return(earliest.possible.days.before.art)
                }
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
                    return(as.numeric(right.boundary - left.boundary))
                } else {
                    return(0)
                }
            }
        )
    }

    return(days.before.art)
}


graft.acute.phase <- function(
    vl.data,
    acute.phase,
    acute.setpoint=1000
) {
    subjects <- names(vl.data)
    regimes <- names(acute.phase)
    grafted.vls <- list()
    for (subject in subjects) {
        vl <- vl.data[[subject]]$vl
        infection.date <- vl.data[[subject]]$infection.date

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


solve.odes <- function(vl.data, grafted.vls) {
    require(deSolve)

    ode.solutions.bin.30 <- list()
    ode.solutions.bin.365 <- list()
    subjects <- names(grafted.vls)
    for (subject in subjects) {
        curr <- vl.data[[subject]]
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
    integration.data,
    ode.solutions.bin.30,
    bin.size=30,  # months, roughly
    possible.half.lives=(1:1800) * 30,
    mc.cores=3
) {
    # Find the decay rate that maximizes the likelihood for each individual.
    require(parallel)
    all.log.likelihoods <- list()  # this will be triply-indexed by pid, regime, and collection date

    subjects <- names(integration.data)
    for (subject in subjects) {
        all.log.likelihoods[[subject]] <- list()
        curr.data <- integration.data[[subject]]

        regimes <- names(ode.solutions.bin.30[[subject]])
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


compute.mles <- function(
    all.log.likelihoods,
    ode.solutions.bin.30,
    bin.size,
    possible.half.lives
) {
    mles <- NULL
    subjects <- names(all.log.likelihoods)
    for (subject in subjects) {
        regimes <- names(all.log.likelihoods[[subject]])
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
        }
    }
    return(mles)
}


compute.bayes.factors <- function(
    integration.data,
    ode.solutions.bin.30,
    all.log.likelihoods
) {
    bayes.factors <- NULL
    subjects <- names(integration.data)
    for (subject in subjects) {
        regimes <- names(ode.solutions.bin.30[[subject]])
        for (regime in regimes) {
            curr.data <- integration.data[[subject]]
            reservoir.dist <- ode.solutions$bin.30[[subject]][[regime]]

            no.decay.ll.no.factorial <- log.likelihood.no.factorial(
                curr.data$days.before.art,
                reservoir.dist$bin.dist.no.decay,
                30
            )
            no.decay.ll <- no.decay.ll.no.factorial + sum(log(1:nrow(curr.data)))

            lls <- all.log.likelihoods[[subject]][[regime]]
            max.idx <- which.max(lls)
            max.ll <- lls[max.idx]

            bayes.factor <- exp(max.ll - no.decay.ll)
            bayes.factors <- rbind(
                bayes.factors,
                data.frame(
                    subject=subject,
                    regime=regime,
                    bayes.factor=bayes.factor
                )
            )
        }
    }
    return(bayes.factors)
}


# We use this helper to compute the breakpoints for making the histograms of
# "days before ART".  This makes sure that the breakpoints cover all of our 
# values properly.
compute.bin.breakpoints <- function(
    days.before.art,  # a vector of integration dates, usually, in "days before ART"
    days.pre.therapy,  # total number of days untreated
    bin.size=365
) {
    breakpoints = seq(0, max(c(days.before.art, days.pre.therapy)), by=bin.size)
    if (!(max(days.before.art) %in% breakpoints)) {
        breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + bin.size)
    }
    return(breakpoints)
}


# This helper performs the GLM regression needed for our model-free estimates.
compute.model.free.estimate <- function(
    days.before.art,
    infection.date,
    art.initiation,
    bin.size=365
) {
    days.pre.therapy <- 
        as.numeric(
            art.initiation - infection.date,
            units="days"
        )

    breakpoints <- compute.bin.breakpoints(
        days.before.art,
        days.pre.therapy,
        bin.size=bin.size
    )
    actual.freqs <- hist(
        days.before.art,
        breaks=breakpoints,
        plot=FALSE
    )
    regression.frame <- data.frame(
        x=1:length(actual.freqs$counts),
        y=actual.freqs$counts
    )
    decay.rate.regression <- glm(
        y ~ x,
        family=poisson,
        data=regression.frame
    )

    return(
        list(
            decay.rate.regression=decay.rate.regression,
            actual.freqs=actual.freqs
        )
    )
}


model.free.half.life <- function(decay.rate.regression) {
    decay.rate.summary <- summary(decay.rate.regression)
    x.coef <- decay.rate.summary$coefficients[2, 1]
    x.se <- decay.rate.summary$coefficients[2, 2]

    half.life <- Inf
    if (x.coef < 0) {
        half.life <- - log(2) / x.coef
    }

    lower.bound <- Inf
    lower.bound.denom <- x.coef - 1.96 * x.se
    if (lower.bound.denom < 0) {
        lower.bound <- - log(2) / lower.bound.denom
    }

    upper.bound <- Inf
    upper.bound.denom <- x.coef + 1.96 * x.se
    if (upper.bound.denom < 0) {
        upper.bound <- - log(2) / upper.bound.denom
    }
    
    return(
        list(
            half.life=half.life,
            lower.bound=lower.bound,
            upper.bound=upper.bound
        )
    )
}
