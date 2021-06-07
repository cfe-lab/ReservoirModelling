# This script produces a plot of the "typical" pVL case for time spans of 3 and 7 years.

# Typically you will have already run the "ode" part.
# load("typical_pvl_demonstration_2021_04_07.RData")

source("reservoir_helpers.r")

# Produce a plot of the viral load curve with three lines:
# - the full viral load curve, uninterrupted
# - the viral load curve, dropped to undetectable at 3 years
# - the viral load curve, dropped to undetectable at 7 years

x.values <- ode.bin.365$solution$time[1:3650000]

pdf("vl_curve_3yr.pdf")

plot(
    x.values, 
    1000 * ode.bin.365$solution$V[1:3650000],
    xlab="days post-infection",
    ylab="Viral load",
    main="Modelled progression of viral load",
    type="l",
    log="y",
    lwd=10,
    col="grey"
)

colour.3.years <- rgb(0, 0, 1, 0.5)
lines(
    x.values,
    c(1000 * ode.bin.365$solution$V[1:(365000 * 3)], rep(30, 365000 * 7)),
    type="l",
    lwd=5,
    lend="butt",
    col=colour.3.years
)

abline(v=x.values[365000 * 3], lty="dashed", lwd=3)

text(
    x=x.values[365000 * 3],
    y=1000000,
    labels=paste("3 years"),
    pos=4
)

legend(
    "topright",
    legend=c("No ART", "ART at 3 years"),
    fill=c("grey", colour.3.years)
)

dev.off()


pdf("vl_curve_7yr.pdf")

plot(
    x.values, 
    1000 * ode.bin.365$solution$V[1:3650000],
    xlab="days post-infection",
    ylab="Viral load",
    main="Modelled progression of viral load",
    type="l",
    log="y",
    lwd=10,
    col="grey"
)

colour.7.years <- rgb(1, 0, 1, 0.5)
lines(
    x.values,
    c(1000 * ode.bin.365$solution$V[1:(365000 * 7)], rep(30, 365000 * 3)),
    type="l",
    lwd=5,
    lend="butt",
    col=colour.7.years
)

abline(v=x.values[365000 * 7], lty="dashed", lwd=3)

text(
    x=x.values[365000 * 7],
    y=1000000,
    labels=paste("7 years"),
    pos=4
)

legend(
    "topright",
    legend=c("No ART", "ART at 7 years"),
    fill=c("grey", colour.7.years)
)

dev.off()


# For each case, plot the hypothetical latent reservoir compositions for
# - no decay
# - 44mo decay
# - 140mo decay
# - 180day decay

for (num.untreated.years in c(3, 7)) {
    bin.freqs <- ode.bin.365$bin.freqs[1:num.untreated.years]

    undecayed <- bin.freqs / sum(bin.freqs)
    dist.44mo.decay <- decay.distribution(bin.freqs, 44 * 30, 365)
    dist.140mo.decay <- decay.distribution(bin.freqs, 140 * 30, 365)
    dist.180day.decay <- decay.distribution(bin.freqs, 180, 365)

    breakpoints = seq(0, 365 * num.untreated.years, by=365)

    max.y <- max(
        c(
            undecayed,
            dist.44mo.decay$bin.dist, 
            dist.140mo.decay$bin.dist, 
            dist.180day.decay$bin.dist 
        )
    )

    pdf(paste0("composition_", num.untreated.years, "yr_untreated.pdf"))

    plot(
        c(0, length(bin.freqs)),
        c(0, max.y),
        main=paste(
            "Reservoir composition,",
            num.untreated.years,
            "years untreated"
        ),
        xlab="Year prior to ART initiation",
        ylab="Proportion",
        type="n",
        xaxt="n"
    )

    axis(
        1,  # this is the x axis
        at=seq(1, length(bin.freqs)) - 0.5,
        labels=seq(length(bin.freqs), 1, by=-1)
    )

    x.coords <- c(0, seq(length(bin.freqs) - length(dist.140mo.decay$bin.dist) + 1, length(bin.freqs)))

    lines(
        x.coords,
        c(undecayed[1], undecayed),
        type="S",
        col="purple",
        lwd=3
    )

    lines(
        x.coords,
        c(dist.140mo.decay$bin.dist[1], dist.140mo.decay$bin.dist),
        type="S",
        col="grey",
        lwd=3
    )

    lines(
        x.coords,
        c(dist.44mo.decay$bin.dist[1], dist.44mo.decay$bin.dist),
        type="S",
        col="black",
        lty="dotdash",
        lwd=3
    )

    lines(
        x.coords,
        c(dist.180day.decay$bin.dist[1], dist.180day.decay$bin.dist),
        type="S",
        col="green",
        lwd=3
    )

    # Build up the legend, with some customization for the "combined" cases.
    legend.labels <- c(
        "no decay",
        "140mo decay",
        "44mo decay",
        "180 day decay"
    )
    legend.colours <- c(
        "purple",
        "grey",
        "black",
        "green"
    )
    legend.ltys <- c(
        "solid",
        "solid",
        "dotdash",
        "solid"
    )
    legend.lwds <- c(3, 3, 3, 3)

    legend(
        "top",
        legend=legend.labels,
        col=legend.colours,
        lty=legend.ltys,
        lwd=legend.lwds
    )
    dev.off()
}


# Next, plot each individual's "actual" pVL curve against the modelled curve.

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

for (pid in unique(vl.data$pid)) {
    pid.vl.data <- vl.data[vl.data$pid == pid,]

    infection.date <- pid.vl.data$est.infection.date[1]
    art.initiation <- pid.vl.data$art.start.date[1]
    max.time <- as.numeric(art.initiation - infection.date, units="days")

    actual.vl.x <- c(0, 1, pid.vl.data$days.after.infection)
    actual.vl.y <- c(0, 30, pid.vl.data$vl)

    typical.vl.x <- ode.bin.365$solution$time[1:(1000 * max.time)]
    typical.vl.y <- 1000 * ode.bin.365$solution$V[1:(1000 * max.time)]

    max.y <- max(actual.vl.y, typical.vl.y)

    pdf(paste0("vl_", pid, ".pdf"))

    plot(
        c(0, max.time),
        c(1, max.y),
        xlab="days post-infection",
        ylab="Viral load",
        main=paste0("Viral load (", pid, ")"),
        type="n",
        log="y"
    )

    lines(
        typical.vl.x,
        typical.vl.y,
        lwd=10,
        col="grey"
    )

    lines(
        actual.vl.x,
        actual.vl.y,
        lwd=3,
        col="red"
    )

    legend(
        "bottomright",
        legend=c("typical", "actual"),
        col=c("grey", "red"),
        lwd=c(10, 3)
    )

    dev.off()
}
