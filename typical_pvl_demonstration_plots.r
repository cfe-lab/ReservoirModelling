# Produce plots of the typical pVL curve, as well as plots of
# the other components of the dynamical system.

# Typically you will have already run the "ode" part.
# load("typical_pvl_demonstration_ode.RData")

source("reservoir_helpers.r")

# Produce a plot of the viral load curve with three lines:
# - the full viral load curve, uninterrupted
# - the viral load curve, dropped to undetectable at 3 years
# - the viral load curve, dropped to undetectable at 7 years

# To keep the plot from being too slow, let's slim down the number of points.
# We want the most detail in the half-year or so and then we can ease up.
time.indices <- c(1:200000, seq(from=201000, 3650000, by=1000))
x.values <- ode.bin.365$solution$time[time.indices]
art.1.colour <- "#5581B0"

for (untreated.years in c(3, 7)) {
    pdf(paste0("vl_curve_", untreated.years, "yr.pdf"))

    par(mar=c(6.5, 6.5, 2, 2) + 0.1)
    plot(
        x.values, 
        1000 * ode.bin.365$solution$V[time.indices],
        xlab=NA,
        ylab=NA,
        xaxt="n",
        yaxt="n",
        type="l",
        log="y",
        lwd=10,
        col="grey",
        cex.lab=3,
        cex.axis=2
    )

    title(
        xlab="Years post-infection",
        line=3.5,
        cex.lab=2.25
    )

    title(
        ylab="Viral load (log 10)",
        line=3.5,
        cex.lab=3
    )

    axis(
        1,  # this is the x axis
        at=seq(0, 365 * 10, by=365),
        # labels=c(0, NA, NA, 3, NA, NA, 6, NA, NA, 9, NA),
        labels=0:10,
        cex.axis=2
    )

    axis(
        2,  # this is the y axis
        at=c(100, 10000, 1000000),
        labels=c(2, 4, 6),
        cex.axis=2
    )

    treated.colour <- rgb(0, 0, 1, 0.5)
    time.indices.before <- time.indices[time.indices <= 365000 * untreated.years]
    time.indices.after <- time.indices[time.indices > 365000 * untreated.years]
    lines(
        x.values,
        c(1000 * ode.bin.365$solution$V[time.indices.before], rep(30, length(time.indices.after))),
        type="l",
        lwd=5,
        lend="butt",
        col=treated.colour
    )

    abline(v=365 * untreated.years, lty="dashed", lwd=3)

    # Add a point representing a hypothetical sampling time one year after this.
    points(
        x=365 * (untreated.years + 1),
        y=60,
        pch=25,
        cex=3,
        lwd=3,
        col=art.1.colour,
        bg=art.1.colour
    )

    text(
        x=365 * untreated.years,
        y=100000,
        labels=paste(untreated.years, "years"),
        pos=4,
        cex=2
    )

    legend(
        "topright",
        legend=c(
            "No ART", 
            paste("ART at", untreated.years, "years"),
            "Proviral sampling"
        ),
        col=c("grey", treated.colour, art.1.colour),
        lty=c("solid", "solid", NA),
        lwd=c(10, 5, 3),
        pch=c(NA, NA, 25),
        pt.cex=c(1, 1, 3),
        pt.bg=c("grey", treated.colour, art.1.colour),
        cex=1.5,
        bg="white"
    )

    dev.off()
}


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

    max.y <- 1

    pdf(paste0("composition_", num.untreated.years, "yr_untreated.pdf"))
    par(mar=c(6.5, 6.5, 2, 2) + 0.1)

    plot(
        c(0, length(bin.freqs)),
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
        at=seq(1, length(bin.freqs)) - 0.5,
        labels=seq(length(bin.freqs), 1, by=-1),
        cex.axis=2
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
        col="red",
        lwd=3
    )

    # Build up the legend, with some customization for the "combined" cases.
    legend.labels <- c(
        "created by viral seeding: no decay",
        "with 140 mo decay",
        "with 44 mo decay",
        "with 6 mo decay"
    )
    legend.colours <- c(
        "purple",
        "grey",
        "black",
        "red"
    )
    legend.ltys <- c(
        "solid",
        "solid",
        "dotdash",
        "solid"
    )
    legend.lwds <- c(3, 3, 3, 3)

    legend(
        "topleft",
        legend=legend.labels,
        col=legend.colours,
        lty=legend.ltys,
        lwd=legend.lwds,
        bty="n",
        cex=1.5
    )
    dev.off()
}
