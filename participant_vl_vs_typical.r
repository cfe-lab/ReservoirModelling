
# Plot each individual's "actual" pVL curve against the modelled curve.

# Typically you will have already run the "ode" part.
# load("typical_pvl_demonstration_ode.RData")

source("reservoir_helpers.r")
source("ode_based_helpers.r")


vl.data <- load.vl.data()


# To keep the plot from being too slow, let's slim down the number of points.
# We want the most detail in the half-year or so and then we can ease up.
time.indices <- c(1:200000, seq(from=201000, 3650000, by=1000))

for (pid in unique(vl.data$pid)) {
    pid.vl.data <- vl.data[vl.data$pid == pid,]

    infection.date <- pid.vl.data$est.infection.date[1]
    art.initiation <- pid.vl.data$art.start.date[1]
    max.time <- as.numeric(art.initiation - infection.date, units="days")

    actual.vl.x <- c(0, 1, pid.vl.data$days.after.infection)
    actual.vl.y <- c(0, 30, pid.vl.data$vl)

    time.indices.pre.art <- time.indices[time.indices < 1000 * max.time]

    typical.vl.x <- ode.bin.365$solution$time[time.indices.pre.art]
    typical.vl.y <- 1000 * ode.bin.365$solution$V[time.indices.pre.art]

    max.y <- 100000000

    pdf(paste0("vl_", pid, ".pdf"))
    par(mar=c(6.5, 6.5, 2, 2) + 0.1)

    plot(
        c(0, max.time),
        c(1, max.y),
        xlab=NA,
        ylab=NA,
        # main=paste0("Viral load (", pid, ")"),
        type="n",
        log="y",
        yaxt="n",
        cex.axis=2
    )

    title(
        xlab="Days post-infection",
        line=3.5,
        cex.lab=2.25
    )

    title(
        ylab="Viral load (log 10)",
        line=3.5,
        cex.lab=3
    )

    axis(
        2,  # this is the y axis
        at=c(100, 10000, 1000000),
        labels=c(2, 4, 6),
        cex.axis=2
    )

    lines(
        typical.vl.x,
        typical.vl.y,
        lwd=10,
        col="grey"
    )

    # Add a special symbol for any data points that are estimated.
    estimated.idx <- which(pid.vl.data$is.estimated.peak != "N")
    points(
        x=pid.vl.data$days.after.infection[estimated.idx],
        y=pid.vl.data$vl[estimated.idx],
        col="purple",
        lwd=3,
        cex=3
    )

    lines(
        actual.vl.x,
        actual.vl.y,
        type="b",
        lty="dashed",
        lwd=3,
        col="red"
    )

    legend(
        "bottomright",
        legend=c("typical", "actual", "estimated peak"),
        col=c("grey", "red", "purple"),
        lty=c("solid", "dashed", NA),
        pch=c(NA, 1, 1),
        lwd=c(10, 3, 3),
        pt.cex=c(1, 1, 3),
        bty="n",
        cex=1.5
    )

    text(
        x=max.time / 2,
        y=50000000,
        labels=pid,
        cex=2
    )

    dev.off()
}
