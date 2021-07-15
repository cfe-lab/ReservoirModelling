source("reservoir_helpers.r")


load.vl.data <- function(
    file.name="../data/BrooksDataVLs.csv"
) {
    vl.data <- read.csv(file.name)
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
    
    return(vl.data)
}


#### 
# Helpers that load integration data.

# This will take a data frame with columns
# - art.start.date
# - est.infection.date
# - integration.date.est
# and return a data frame with columns
# - untreated.period
# - days.before.art.raw (number of days prior to ART initiation that each sequence is estimated
#   to have been introduced to the reservoir)
# - days.before.art (the above, but negative values are truncated to 0)
compute.day.columns <- function(integration.data) {
    untreated.period <- as.numeric(
        integration.data$art.start.date - integration.data$est.infection.date,
        units="days"
    )
    days.before.art.raw <- as.numeric(
        integration.data$art.start.date - integration.data$integration.date.est,
        units="days"
    )
    days.before.art <- sapply(
        days.before.art.raw,
        function (x) {
            if (x >= 0) {
                return(x)
            }
            return(0)
        }
    )
    return(
        data.frame(
            untreated.period=untreated.period,
            days.before.art.raw=days.before.art.raw,
            days.before.art=days.before.art
        )
    )
}


integration.data.column.names <- c(
    "pid",
    "est.infection.date",
    "art.start.date",
    "seq.id",
    "collection.date",
    "integration.date.est",
    "integration.date.lower",
    "integration.date.upper"
)
load.integration.data <- function(
    file.name="../data/BrooksData.csv"
) {
    brooks.data <- read.csv(file.name)[, 1:8]
    names(brooks.data) <- integration.data.column.names
    for (col.idx in c(2, 3, 5, 6, 7, 8)) {
        brooks.data[[col.idx]] <- strptime(brooks.data[[col.idx]], "%d-%b-%y")
    }

    return(
        cbind(
            brooks.data,
            compute.day.columns(brooks.data)
        )
    )
}


# This function will load the data *with* the duplicates.
load.integration.data.with.duplicates <- function(
    with.duplicates.file.name="../data/BrooksDataWithDuplicates.csv"
) {
    with.duplicates.data <- read.csv(with.duplicates.file.name)[, 1:6]
    names(with.duplicates.data) <- integration.data.column.names[1:6]
    for (col.idx in c(2, 3, 5, 6)) {
        with.duplicates.data[[col.idx]] <- strptime(with.duplicates.data[[col.idx]], "%d-%b-%y")
    }

    return(
        cbind(
            with.duplicates.data,
            compute.day.columns(with.duplicates.data)
        )
    )
}


# This will combine the results from the previous two functions.
load.combined.integration.data <- function(without, with) {
    useful.columns <- c(1:6, 9:11)

    pids.with.duplicates <- unique(with$pid)

    combined.df <- NULL
    for (pid in unique(without$pid)) {
        curr.pid.df <- without[without$pid == pid, useful.columns]
        if (pid %in% pids.with.duplicates) {
            curr.pid.df <- with[with$pid == pid,]
        }
        combined.df <- rbind(combined.df, curr.pid.df)
    }

    return(combined.df)
}


date.sanity.check <- function(vl.data, sample.data) {
    for (pid in unique(vl.data$pid)) {
        pid.vl.data <- vl.data[vl.data$pid == pid,]

        infection.date <- pid.vl.data$est.infection.date[1]
        art.initiation <- pid.vl.data$art.start.date[1]

        pid.sample.data <- sample.data[sample.data$pid == pid,]
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
    }
}


perform.odes <- function(vl.data) {
    ode.known.bin.30 <- list()
    ode.known.bin.365 <- list()
    ode.typical.bin.30 <- list()
    ode.typical.bin.365 <- list()
    for (pid in unique(vl.data$pid)) {
        pid.vl.data <- vl.data[vl.data$pid == pid,]

        infection.date <- pid.vl.data$est.infection.date[1]
        art.initiation <- pid.vl.data$art.start.date[1]

        # Compute the ODE for this individual, and bin the reservoir into 30 day intervals.
        cat("Solving ODEs for ", pid, ".\n", sep="")
        ode.known.bin.30[[pid]] <- undecayed.reservoir.distribution.given.vl(
            c(0, pid.vl.data$days.after.infection),
            c(0, pid.vl.data$vl / 1000),  # the VL was given in copies/mL; convert to copies/uL
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

    return(
        list(
            bin.30=ode.solutions.bin.30,
            bin.365=ode.solutions.bin.365
        )
    )
}


# This function will perform the likelihood calculations
# and compute the MLEs.  The result will be a list with two entries:
#
# - all.log.likelihoods
# This is a list that will be triply-indexed by:
# VL info: "known" or "typical"
# PID
# collection date (this could be "combined" also)
#
# - mles
# This will be a dataframe with columns
# vl.info ("known" or "typical")
# pid
# col.date (for most, there's only one collection date, but for those with two collection
# dates, there will be separate entries for both, as well as one for "combined")
# mle
# lower.bound
# upper.bound
ode.based.likelihood.estimates <- function(
    integration.data,
    ode.solutions.bin.30,
    bin.size=30,
    possible.half.lives=(1:200) * 30,
    parallel.cores=3
) {
    require(parallel)

    all.log.likelihoods <- list()
    mles <- NULL

    for (vl.info in c("known", "typical")) {
        all.log.likelihoods[[vl.info]] <- list()

        for (curr.pid in unique(integration.data$pid)) {
            cat(
                "Processing data for ",
                curr.pid,
                "....\n",
                sep=""
            )
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
                lls <- (
                    unlist(log.likelihoods.no.factorial) + 
                    sum(log(1:nrow(curr.data)))
                )
                all.log.likelihoods[[vl.info]][[curr.pid]][[col.date]] <- lls

                # Record the MLE and error bounds, calculated using the Fisher 
                # information matrix.
                max.idx <- which.max(lls)
                mle <- possible.half.lives[max.idx]
                estimated.variance <- (
                    1 / fisher.information(
                        max.idx, 
                        ode.solutions.bin.30[[vl.info]][[curr.pid]]$bin.freqs
                    )
                )
                lower.bound <- max(0, mle - 2 * sqrt(estimated.variance) * bin.size)
                upper.bound <- mle + 2 * sqrt(estimated.variance) * bin.size
                mles <- rbind(
                    mles,
                    data.frame(
                        vl.info=vl.info,
                        pid=curr.pid,
                        col.date=col.date,
                        mle=mle,
                        lower.bound=lower.bound,
                        upper.bound=upper.bound
                    )
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

    return(
        list(
            all.log.likelihoods=all.log.likelihoods,
            mles=mles
        )
    )
}


plot.log.likelihoods <- function(
    lls,
    possible.half.lives,
    plot.title,
    mle,
    lower.bound,
    upper.bound
) {
    plot(
        possible.half.lives,
        lls,
        main=plot.title,
        xlab="Reservoir half life (days)",
        ylab="Log likelihood"
    )

    # Plot the MLE and its error bounds.
    abline(v=mle, lty="dashed")
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
}


# The colour scheme:
# seroconversion     #BA2A19
# 1 Year             #000000
# Last ART-naive     #1332F5
# ART 1              #5581B0
# ART 2              #74DCD0
default.art.1.colour <- "#5581B0"
default.art.2.colour <- "#74DCD0"
plot.ode.based.reservoir.composition <- function(
    integration.days.before.art,  # integration dates, in "days before ART"
    collection.dates,  # actual dates that represent the dates the proviruses were sampled
    reservoir.dist.bin.freqs,  # undecayed latent reservoir age distribution, binned by year
    days.pre.therapy,  # how long the person's untreated period was
    mle,
    lower.bound,
    upper.bound,
    show.legend=FALSE,
    art.1.colour=default.art.1.colour,
    art.2.colour=default.art.2.colour
) {
    dist.44mo.decay <- decay.distribution(reservoir.dist.bin.freqs, 44 * 30, 365)
    dist.140mo.decay <- decay.distribution(reservoir.dist.bin.freqs, 140 * 30, 365)
    dist.best.fit <- decay.distribution(reservoir.dist.bin.freqs, mle, 365)

    breakpoints = seq(0, max(c(integration.days.before.art, days.pre.therapy)), by=365)
    if (!(max(integration.days.before.art) %in% breakpoints)) {
        breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + 365)
    }
    actual.freqs <- hist(
        integration.days.before.art,
        breaks=breakpoints,
        plot=FALSE
    )
    emp.dist <- actual.freqs$counts / sum(actual.freqs$counts)

    max.y <- 1

    par(mar=c(6.5, 6.5, 2, 2) + 0.1)
    plot(
        c(0, length(emp.dist)),
        c(0, max.y),
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
        ylab="Proportion of proviruses",
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
            col=art.1.colour,
            lend="butt"
        )
    } else {
        # Plot the empirical distribution as a stacked bar graph.  
        # There are at most two collection dates.
        distinct.dates <- unique(collection.dates)
        first.collected <- integration.days.before.art[collection.dates == min(distinct.dates)]
        last.collected <- integration.days.before.art[collection.dates == max(distinct.dates)]

        first.freqs <- hist(
            first.collected,
            breaks=breakpoints,
            plot=FALSE
        )
        last.freqs <- hist(
            last.collected,
            breaks=breakpoints,
            plot=FALSE
        )

        total <- sum(first.freqs$counts, last.freqs$counts)

        # First plot the earlier ones as with the other plots.
        lines(
            seq(length(emp.dist), 1, by=-1) - 0.5,
            first.freqs$counts / total,
            type="h",
            lwd=20,
            col=art.1.colour,
            lend="butt"
        )

        # Then plot the later ones on top.
        segments(
            x0=seq(length(emp.dist), 1, by=-1) - 0.5,
            y0=first.freqs$counts / total,
            y1=(first.freqs$counts + last.freqs$counts) / total,
            lwd=20,
            col=art.2.colour,
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
        col="red",
        lwd=3
    )

    # Whether there is a legend or not determines where we plot the half-life 
    # estimate and the 95% CI.
    est.text.x <- 0
    est.text.y <- max.y - 0.015

    if (show.legend) {
        # Build up the legend, with some customization for the "combined" cases.                
        legend.labels <- "sampled"
        legend.colours <- art.1.colour
        legend.ltys <- "solid"
        legend.lwds <- 20
        if (col.date == "combined") {
            legend.labels <- c(
                "sampling 1",
                "sampling 2"
            )
            legend.colours <- c(art.1.colour, art.2.colour)
            legend.ltys <- c("solid", "solid")
            legend.lwds <- c(20, 20)
        }
        legend.labels <- c(
            legend.labels,
            "140 mo decay",
            "44 mo decay",
            "best-fit decay"
        )
        legend.colours <- c(
            legend.colours,
            "grey",
            "black",
            "red"
        )
        legend.ltys <- c(
            legend.ltys,
            "solid",
            "dotdash",
            "solid"
        )
        legend.lwds <- c(legend.lwds, 3, 3, 3)

        legend.location <- legend(
            "topleft",
            legend=legend.labels,
            col=legend.colours,
            lty=legend.ltys,
            lwd=legend.lwds,
            bty="n",
            cex=1.5
        )
        est.text.x <- legend.location$rect$left
        est.text.y <- legend.location$rect$top - legend.location$rect$h - 0.015
        est.text <- substitute(
            paste(
                t[1/2],
                " = ",
                mle.formatted,
                " yr, 95% CI [",
                lb.formatted,
                ", ",
                ub.formatted,
                ")",
                sep=""
            ),
            list(
                mle.formatted=round(mle / 365, digits=2),
                lb.formatted=round(lower.bound / 365, digits=2),
                ub.formatted=round(upper.bound / 365, digits=2)
            )
        )
        text(
            x=est.text.x,
            y=est.text.y,
            labels=est.text,
            pos=4,
            cex=1.5,
            offset=0.5
        )
    } else {
        est.text <- substitute(
            paste(
                t[1/2],
                " = ",
                mle.formatted,
                " yr",
                sep=""
            ),
            list(mle.formatted=round(mle / 365, digits=2))
        )

        text(
            x=est.text.x,
            y=est.text.y,
            labels=est.text,
            pos=4,
            cex=1.5,
            offset=0.5
        )

        ci.text <- substitute(
            paste(
                "95% CI [",
                lb.formatted,
                ", ",
                ub.formatted,
                ")",
                sep=""
            ),
            list(
                lb.formatted=round(lower.bound / 365, digits=2),
                ub.formatted=round(upper.bound / 365, digits=2)
            )
        )

        text(
            x=est.text.x,
            y=est.text.y - 0.062,
            labels=ci.text,
            pos=4,
            cex=1.5,
            offset=0.5
        )
    }

    text(
        x=length(emp.dist) * 0.8,
        y=0.95,
        labels=pid,
        cex=2
    )
}
