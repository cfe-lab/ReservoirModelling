# Follow the recipe in ode.r, but with a known viral load curve.

library(deSolve)

# These constants are taken from the literature.
days.to.peak <- 23
peak.vl <- 12902  # copies/mL
days.to.undetectable <- 154

# Read in the data.
all.subjects <- list()

for (subject in c("p1", "p2", "p3", "p4")) {
    # The VL CSV file has columns (date, viral.load, comments).
    vl.csv <- read.csv(paste(subject, "/VL-VL.csv", sep=""))
    vl.csv$date <- strptime(vl.csv$date, "%d/%m/%Y")

    peak.csv <- read.csv(paste(subject, "/Peak-Peak.csv", sep=""))
    infection.date <- strptime(peak.csv$est.infection.date[1], "%d/%m/%Y")
    art.initiation <- strptime(peak.csv$art.initiation[1], "%d/%m/%Y")

    all.subjects[[subject]] <- list(
        vl=vl.csv,
        infection.date=infection.date,
        art.initiation=art.initiation
    )
}

acute.time <- c(
    0,
    days.to.peak,
    days.to.peak + days.to.undetectable
)
acute.vls <- c(0, peak.vl, 20)
acute.comments <- rep("graft", 3)
graft.df <- data.frame(
    time=acute.time,
    viral.load=acute.vls,
    comments=acute.comments
)
grafted.vls <- list()
for (subject in c("p1", "p2", "p3", "p4")) {
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

    grafted.vl <- NULL
    if (acute.time[3] > time[1]) {
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
    grafted.vls[[subject]] <- grafted.vl
}


# Next, we define the system of ODEs.
source("reservoir_helpers.r")

state.vl.known <- state[1:6]

all.odes <- list()
for (subject in c("p1", "p2", "p3", "p4")) {
    vl <- grafted.vls[[subject]]
    virus <- approxfun(vl$time, vl$viral.load, rule=2)

    derivatives.vl.known <- function(t, state, parameters) {
        with(
            as.list(c(state, parameters)), 
            {
                V <- virus(t)

                dS <- alpha.S - delta.S * S - beta * S * V
                dA.P <- (1 - lambda) * tau * beta * S * V - delta.I * A.P - kappa * A.P * E
                dA.U <- (1 - lambda) * (1 - tau) * beta * S * V - delta.I * A.U - kappa * A.U * E
                dL.P <- lambda * tau * beta * S * V
                dL.U <- lambda * (1 - tau) * beta * S * V
                dE <- alpha.E + omega * (A.P + A.U) * E / (E + E.50) - delta.E * E

                return(
                    list(
                        c(dS, dA.P, dA.U, dL.P, dL.U, dE),
                        V=V
                    )
                )
            }
        )
    }

    infection.date <- all.subjects[[subject]]$infection.date
    art.initiation <- all.subjects[[subject]]$art.initiation
    max.time <- as.numeric(art.initiation - infection.date, units="days")
    times <- seq(0, max.time, by=0.001)

    cat("Solving ODEs (VL known) for subject ", subject, "....\n", sep="")
    run.time <- system.time(
        vl.known <- ode(
            y=state.vl.known,
            times=times,
            func=derivatives.vl.known,
            parms=parameters,
            maxsteps=25000
        )
    )
    cat("Time elapsed:\n")
    print(run.time)

    cat("Solving ODEs (idealized) for subject ", subject, "....\n", sep="")
    run.time <- system.time(
        idealized <- ode(
            y=state,
            times=times,
            func=derivatives,
            parms=parameters,
            maxsteps=25000
        )
    )
    cat("Time elapsed:\n")
    print(run.time)

    all.odes[[subject]] <- list(vl.known=vl.known, idealized=idealized)
}


# For each subject, make a series of plots for a few different decay rates comparing
# both the idealized and known VL results.
decay.rates <- c(14, 180, 365, 44 * 30, 140 * 30)
decay.rate.strings <- c("2wk", "6mo", "1yr", "44mo", "140mo")
for (subject in c("p1", "p2", "p3", "p4")) {
    idealized.df <- as.data.frame(all.odes[[subject]]$idealized)
    vl.known.df <- as.data.frame(all.odes[[subject]]$vl.known)

    idealized.bin.freqs <- bin.helper(idealized.df, bin.size=365, grid.size=0.001)
    vl.known.bin.freqs <- bin.helper(vl.known.df, bin.size=365, grid.size=0.001)

    x.labels <- length(idealized.bin.freqs) - 1:length(idealized.bin.freqs) + 1
    
    barplot.matrix <- rbind(
        idealized.bin.freqs / sum(idealized.bin.freqs), 
        vl.known.bin.freqs / sum(vl.known.bin.freqs)
    )
    pdf(paste(subject, "_no_decay.pdf", sep=""))
    barplot(
        barplot.matrix,
        names.arg=x.labels,
        main=paste("Composition of", subject, "reservoir by age (no decay)"),
        xlab="Year prior to ART initiation",
        ylab="Proportion of reservoir",
        beside=TRUE,
        legend.text=c("idealized", "VL known")
    )
    dev.off()
        
    for (idx in 1:length(decay.rates)) {
        decay.rate <- decay.rates[idx]
        decay.rate.str <- decay.rate.strings[idx]
        idealized.dist <- decay.distribution(idealized.bin.freqs, decay.rate, bin.size=365)$bin.dist
        vl.known.dist <- decay.distribution(vl.known.bin.freqs, decay.rate, bin.size=365)$bin.dist

        barplot.matrix <- rbind(idealized.dist, vl.known.dist)
        pdf(paste(subject, "_", decay.rate.str, "_decay.pdf", sep=""))
        barplot(
            barplot.matrix,
            names.arg=x.labels,
            main=paste("Composition of ", subject, " reservoir by age (", decay.rate.str, " decay)", sep=""),
            xlab="Year prior to ART initiation",
            ylab="Proportion of reservoir",
            beside=TRUE,
            legend.text=c("idealized", "VL known")
        )
        dev.off()
    }
}
