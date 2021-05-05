# This script will generate plots analogous to Pankau fig6B, plotting the 
# half-lives we estimated using our model-free method against the QVOA and DNA
# decay rates from the literature.

# The actual estimation will proceed as in our model-free plots, and then we will
# produce a plot showing all of the estimates on the same plot.

# This is commented out as we will often have already loaded it when we run this script.
# load("brooks_data_2021_05_03.RData")


# This largely follows the method of the Apr 13, 2021 model-free script.
vl.info <- "known"
all.regressions <- list()

# Our colour scheme for the "bar" parts of these graphs.
art.1.colour <- "#5581B0"
art.2.colour <- "#74DCD0"

half.lives <- NULL
for (pid in names(all.log.likelihoods[[vl.info]])) {
    all.regressions[[pid]] <- list()

    curr.pid.data <- integration.data[integration.data$pid == pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    for (col.date in names(all.log.likelihoods[[vl.info]][[pid]])) {
        curr.data <- curr.pid.data
        if (col.date != "combined") {
            curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
        }

        breakpoints = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=365)
        if (!(max(curr.data$days.before.art) %in% breakpoints)) {
            breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + 365)
        }

        actual.freqs <- hist(
            curr.data$days.before.art,
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
        all.regressions[[pid]][[col.date]] <- decay.rate.regression

        x.fit.values <- data.frame(x=seq(0.5, length(actual.freqs$counts) + 0.5, by=0.1))
        fit <- predict(decay.rate.regression, newdata=x.fit.values, se.fit=TRUE)

        predictions <- exp(fit$fit)
        upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
        lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)

        # Add text listing the half life.
        decay.rate.summary <- summary(decay.rate.regression)
        x.coef <- decay.rate.summary$coefficients[2, 1]
        x.se <- decay.rate.summary$coefficients[2, 2]
        
        # t_{1/2} = - log(2) / x.coef
        half.life <- - log(2) / x.coef
        half.life.upper <- Inf
        if (x.coef + 1.96 * x.se < 0) {
            half.life.upper <- - log(2) / (x.coef + 1.96 * x.se)
        }
        half.life.lower <- - log(2) / (x.coef - 1.96 * x.se)

        half.lives <- rbind(
            half.lives,
            data.frame(
                pid=pid,
                col.date=col.date,
                half.life=half.life,
                half.life.upper=half.life.upper,
                half.life.lower=half.life.lower
            )
        )
    }
}


# For the first plot, plot all the half-lives alongside the QVOA and DNA half-lives.
# The numerators below are all in months so we divide by 12 to get years.
qvoa.half.life <- 44.2 / 12
qvoa.half.life.upper <- 114.5 / 12
qvoa.half.life.lower <- 27.4 / 12
dna.half.life <- 140.4 / 12
dna.half.life.upper <- Inf
dna.half.life.lower <- 75.6 / 12

# Only include the "combined" half-lives for N133M, Z634F, Z1165M, and Z1788F.
two.sampling.pids <- c("N133M", "Z634F", "Z1165M", "Z1788F")
plot.1.rows <- !(half.lives$pid %in% two.sampling.pids) | (half.lives$col.date == "combined")
plot.1.data <- half.lives[plot.1.rows,]

# Figure out the extents of the plot first.  For entries where the upper bound is infinity,
# we will represent it with an arrow pointing upward.
max.y <- max(
    plot.1.data$half.life.upper[plot.1.data$half.life.upper < Inf], 
    qvoa.half.life.upper, 
    plot.1.data$half.life + 2, 
    dna.half.life + 2
)

pdf("half_lives.pdf")
par(mar=c(8.5, 6.5, 2, 2) + 0.1)
plot(
    x=c(1, nrow(plot.1.data) + 2),
    y=c(0, max.y),
    xlab=NA,
    ylab=NA,
    type="n",
    xaxt="n",
    yaxt="n",
    cex.lab=3,
    cex.axis=2
)

title(
    ylab="Half-life (years)",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=1:(nrow(plot.1.data) + 2),
    labels=c(as.character(plot.1.data$pid), "QVOA", "DNA"),
    cex.axis=2,
    las=2
)
axis(2, cex.axis=2)

# We need to plot the PIDs with two sampling points separately so we can use
# a different symbol (a circle with art.1.colour as the outline and art.2.colour
# as the background).

two.sampling.indices <- plot.1.data$pid %in% two.sampling.pids
points(
    x=which(!two.sampling.indices),
    y=plot.1.data$half.life[!two.sampling.indices],
    cex=3,
    pch=19,
    col=art.1.colour
)
points(
    x=which(two.sampling.indices),
    y=plot.1.data$half.life[two.sampling.indices],
    cex=3,
    lwd=5,
    pch=21,
    col=art.1.colour,
    bg=art.2.colour
)

# Now make the error bars.  We do the top and bottom halves separately.
arrows(
    x0=1:nrow(plot.1.data),
    y0=plot.1.data$half.life,
    y1=plot.1.data$half.life.lower,
    angle=90,
    length=0.1,
    code=2,
    lwd=2
)
upper.bound.inf <- is.infinite(plot.1.data$half.life.upper)
arrows(
    x0=which(!upper.bound.inf),
    y0=plot.1.data$half.life[!upper.bound.inf],
    y1=plot.1.data$half.life.upper[!upper.bound.inf],
    angle=90,
    length=0.1,
    code=2,
    lwd=2
)
arrows(
    x0=which(upper.bound.inf),
    y0=plot.1.data$half.life[upper.bound.inf],
    y1=plot.1.data$half.life[upper.bound.inf] + 2,
    angle=30,
    length=0.25,
    code=2,
    lwd=2
)

points(
    x=c(nrow(plot.1.data) + 1, nrow(plot.1.data) + 2),
    y=c(qvoa.half.life, dna.half.life),
    cex=4,
    pch=18
)
arrows(
    x0=nrow(plot.1.data) + 1,
    y0=qvoa.half.life.lower,
    y1=qvoa.half.life.upper,
    angle=90,
    length=0.1,
    code=3,
    lwd=2
)
arrows(
    x0=nrow(plot.1.data) + 2,
    y0=dna.half.life,
    y1=dna.half.life.lower,
    angle=90,
    length=0.1,
    code=2,
    lwd=2
)
arrows(
    x0=nrow(plot.1.data) + 2,
    y0=dna.half.life,
    y1=dna.half.life + 2,
    angle=30,
    length=0.25,
    code=2,
    lwd=2
)
dev.off()


# Plot 2: plot the half-lives for those individuals with two sampling
# points both separately and with their integration data combined.
plot.2.rows <- half.lives$pid %in% two.sampling.pids
plot.2.data <- half.lives[plot.2.rows,]

first.sampling <- c(1, 4, 7, 10)
second.sampling <- c(2, 5, 8, 11)
combined.sampling <- c(3, 6, 9, 12)

upper.bound.inf <- is.infinite(plot.2.data$half.life.upper)
max.y <- max(
    plot.2.data$half.life.upper[!upper.bound.inf],
    plot.2.data$half.life[upper.bound.inf] + 2
)

pdf("half_lives_two_sampling_points.pdf")
par(mar=c(8.5, 6.5, 2, 2) + 0.1)
plot(
    x=c(1, nrow(plot.2.data)),
    y=c(0, max.y),
    xlab=NA,
    ylab=NA,
    type="n",
    xaxt="n",
    yaxt="n",
    cex.lab=3,
    cex.axis=2
)

title(
    ylab="Half-life (years)",
    line=3.5,
    cex.lab=3
)

axis(
    1,  # this is the x axis
    at=second.sampling,  # only label the middle of each triplet of columns
    labels=as.character(plot.2.data$pid[second.sampling]),
    cex.axis=2,
    las=2
)
axis(2, cex.axis=2)

points(
    x=first.sampling,
    y=plot.2.data$half.life[first.sampling],
    cex=3,
    pch=19,
    col=art.1.colour
)

points(
    x=second.sampling,
    y=plot.2.data$half.life[second.sampling],
    cex=3,
    pch=19,
    col=art.2.colour
)

points(
    x=combined.sampling,
    y=plot.2.data$half.life[combined.sampling],
    cex=3,
    lwd=5,
    pch=21,
    col=art.1.colour,
    bg=art.2.colour
)

arrows(
    x0=which(!upper.bound.inf),
    y0=plot.2.data$half.life.lower[!upper.bound.inf],
    y1=plot.2.data$half.life.upper[!upper.bound.inf],
    angle=90,
    length=0.1,
    code=3,
    lwd=2
)
arrows(
    x0=which(upper.bound.inf),
    y0=plot.2.data$half.life[upper.bound.inf],
    y1=plot.2.data$half.life.lower[upper.bound.inf],
    angle=90,
    length=0.1,
    code=2,
    lwd=2
)
arrows(
    x0=which(upper.bound.inf),
    y0=plot.2.data$half.life[upper.bound.inf],
    y1=plot.2.data$half.life[upper.bound.inf] + 2,
    angle=30,
    length=0.25,
    code=2,
    lwd=2
)


for (separator.loc in c(3.5, 6.5, 9.5)) {
    abline(v=separator.loc, lty="dashed", col="grey")
}

legend(
    "topleft",
    legend=c("sampling 1", "sampling 2", "combined"),
    col=c(art.1.colour, art.2.colour, art.1.colour),
    lty=NA,
    lwd=c(1, 1, 5),
    pch=c(19, 19, 21),
    pt.cex=3,
    pt.bg=c(art.1.colour, art.2.colour, art.2.colour),
    cex=1.5,
    bg="white"
)

dev.off()
