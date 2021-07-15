# Make plots of all of the components of our "typical" dynamical system.

# Typically you will have already run the "ode" part.
# load("typical_pvl_demonstration_ode.RData")

source("reservoir_helpers.r")


# For these plots, we focus on the first 25 weeks.
time.indices <- c(1:(10 * 7 * 1000), seq(70000, 25 * 7 * 1000, by=1000))
x.values <- ode.bin.365$solution$time[time.indices]
x.axis.tick.positions <- seq(0, 7 * 25, by=7 * 5)
x.axis.tick.labels <- seq(0, 25, by=5)
x.axis.cex <- 2


# We'll make plots for:
# S (susceptible cells)
# V (viral load)
# A_P and A_U (these can be plotted together)
# L_P and L_U (these can also be plotted together)
# E (the "adaptive immune component")


pdf("typical_model_vl.pdf")

par(mar=c(6.5, 6.5, 2, 2) + 0.1)
plot(
    x.values, 
    1000 * ode.bin.365$solution$V[time.indices],  # convert cells/uL to cells/mL
    xlab=NA,
    ylab=NA,
    xaxt="n",
    yaxt="n",
    type="l",
    log="y",
    lwd=10,
    col="purple",
    cex.lab=3,
    cex.axis=2,
    cex.main=2,
    main="Viral load"
)

title(
    xlab="Weeks post-infection",
    line=3.5,
    cex.lab=2.25
)

title(
    ylab="cells/mL",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=x.axis.tick.positions,
    labels=x.axis.tick.labels,
    cex.axis=x.axis.cex
)

axis(
    2,  # this is the y axis
    at=c(100, 10000, 1000000),
    labels=c(2, 4, 6),
    cex.axis=2
)

dev.off()


pdf("typical_model_susceptible.pdf")

par(mar=c(6.5, 6.5, 2, 2) + 0.1)
plot(
    x.values, 
    ode.bin.365$solution$S[time.indices],
    xlab=NA,
    ylab=NA,
    xaxt="n",
    type="l",
    lwd=10,
    col="grey",
    cex.lab=3,
    cex.axis=2,
    cex.main=2,
    main="Susceptible"
)

title(
    xlab="Weeks post-infection",
    line=3.5,
    cex.lab=2.25
)

title(
    ylab="cells/uL",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=x.axis.tick.positions,
    labels=x.axis.tick.labels,
    cex.axis=x.axis.cex
)

dev.off()


pdf("typical_model_active.pdf")

par(mar=c(6.5, 6.5, 2, 2) + 0.1)

# This plot covers both A_P and A_U, so we set up the axes first.
max.y <- max(ode.bin.365$solution$A.P, ode.bin.365$solution$A.U)

plot(
    range(x.values),
    range(0, max.y),
    xlab=NA,
    ylab=NA,
    xaxt="n",
    type="n",
    # lwd=10,
    # col="grey",
    cex.lab=3,
    cex.axis=2,
    cex.main=2,
    main="Active"
)

col.p <- col2rgb("pink")
productive.col <- rgb(col.p[1], col.p[2], col.p[3], alpha=127, maxColorValue=255)
col.u <- col2rgb("red")
unproductive.col <- rgb(col.u[1], col.u[2], col.u[3], alpha=127, maxColorValue=255)

lines(
    x.values,
    ode.bin.365$solution$A.P[time.indices],
    xlab=NA,
    ylab=NA,
    xaxt="n",
    yaxt="n",
    lwd=10,
    col=productive.col,
    cex.lab=3,
    cex.axis=2
)
lines(
    x.values,
    ode.bin.365$solution$A.U[time.indices],
    xlab=NA,
    ylab=NA,
    xaxt="n",
    yaxt="n",
    lwd=10,
    col=unproductive.col,
    cex.lab=3,
    cex.axis=2
)

title(
    xlab="Weeks post-infection",
    line=3.5,
    cex.lab=2.25
)

title(
    ylab="cells/uL",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=x.axis.tick.positions,
    labels=x.axis.tick.labels,
    cex.axis=x.axis.cex
)

legend(
    "topright",
    legend=c("productive", "unproductive"),
    col=c(productive.col, unproductive.col),
    lwd=10,
    bty="n",
    cex=1.5
)

dev.off()


pdf("typical_model_latent.pdf")

par(mar=c(6.5, 6.5, 2, 2) + 0.1)

# This plot covers both L_P and L_U, so we set up the axes first.
max.y <- max(ode.bin.365$solution$L.P[time.indices], ode.bin.365$solution$L.U[time.indices])

plot(
    range(x.values),
    range(0, max.y),
    xlab=NA,
    ylab=NA,
    xaxt="n",
    type="n",
    # lwd=10,
    # col="grey",
    cex.lab=3,
    cex.axis=2,
    cex.main=2,
    main="Latent"
)

col.p <- col2rgb("turquoise")
productive.col <- rgb(col.p[1], col.p[2], col.p[3], alpha=127, maxColorValue=255)
col.u <- col2rgb("blue")
unproductive.col <- rgb(col.u[1], col.u[2], col.u[3], alpha=127, maxColorValue=255)

lines(
    x.values,
    ode.bin.365$solution$L.P[time.indices],
    xlab=NA,
    ylab=NA,
    xaxt="n",
    yaxt="n",
    lwd=10,
    col=productive.col,
    cex.lab=3,
    cex.axis=2
)
lines(
    x.values,
    ode.bin.365$solution$L.U[time.indices],
    xlab=NA,
    ylab=NA,
    xaxt="n",
    yaxt="n",
    lwd=10,
    col=unproductive.col,
    cex.lab=3,
    cex.axis=2
)

title(
    xlab="Weeks post-infection",
    line=3.5,
    cex.lab=2.25
)

title(
    ylab="cells/uL",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=x.axis.tick.positions,
    labels=x.axis.tick.labels,
    cex.axis=x.axis.cex
)

legend(
    "right",
    legend=c("productive", "unproductive"),
    col=c(productive.col, unproductive.col),
    lwd=10,
    bty="n",
    cex=1.5
)

dev.off()


pdf("typical_model_adaptive_immune_response.pdf")

par(mar=c(6.5, 6.5, 2, 2) + 0.1)
plot(
    x.values, 
    ode.bin.365$solution$E[time.indices],
    xlab=NA,
    ylab=NA,
    xaxt="n",
    type="l",
    lwd=10,
    col="forestgreen",
    cex.lab=3,
    cex.axis=2,
    cex.main=2,
    main="Adaptive immune response"
)

title(
    xlab="Weeks post-infection",
    line=3.5,
    cex.lab=2.25
)

title(
    ylab="cells/uL",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=x.axis.tick.positions,
    labels=x.axis.tick.labels,
    cex.axis=x.axis.cex
)

dev.off()
