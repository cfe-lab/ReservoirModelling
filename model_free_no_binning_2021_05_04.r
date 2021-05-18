# This performs the model-free fit but with no binning at all on the GLM.
# (More accurately, the dates are binned by *day*, which is all the granularity
# we have in our data anyway.)

# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.

# This is commented out as we will often have already loaded it when we run this script.
# load("brooks_data_2021_05_03.RData")


# This will be doubly-indexed by:
#  - PID
#  - collection date (could also be "combined")
# As the VL data doesn't enter into this analysis, there's nothing to separate the
# "known" and "typical" VL cases.
vl.info <- "known"
all.regressions <- list()
half.lives <- NULL

# Our colour scheme for the "bar" parts of these graphs.
art.1.colour <- "#5581B0"
art.2.colour <- "#74DCD0"
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

        total.count <- sum(actual.freqs$counts)
        # max.y <- max(c(actual.freqs$counts, upper.bounds)) / total.count
        # max.y <- min(1, max.y)
        max.y <- 1

        cairo_pdf(
            paste0(
                "model_free_no_binning_", 
                pid,
                "_",
                col.date,
                ".pdf"
            )
        )
        par(mar=c(6.5, 6.5, 2, 2) + 0.1)
        plot(
            c(0, length(actual.freqs$counts)),
            c(0, max.y),
            xlab=NA,
            ylab=NA,
            type="n",
            xaxt="n",
            yaxt="n",
            cex.lab=3,
            cex.axis=2
        )

        # Some defaults, and then some customization for p3.
        x.label <- "Year prior to ART initiation"
        x.label.cex <- 2.25

        title(
            xlab=x.label,
            line=3.5,
            cex.lab=x.label.cex
        )

        title(
            ylab="Proportion of proviruses",
            line=3.5,
            cex.lab=3
        )

        axis(
            1,  # this is the x axis
            at=seq(1, length(actual.freqs$counts)) - 0.5,
            labels=seq(length(actual.freqs$counts), 1, by=-1),
            cex.axis=2
        )
        axis(2, cex.axis=2)

        if (col.date != "combined") {
            lines(
                seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
                actual.freqs$counts / total.count,
                type="h",
                lwd=20,
                col=art.1.colour,
                lend="butt"
            )
        } else {
            # Plot the empirical distribution as a stacked bar graph.  
            # There are at most two collection dates.
            collection.dates <- unique(curr.data$collection.date)
            first.collected <- curr.data[curr.data$collection.date == min(collection.dates),]
            last.collected <- curr.data[curr.data$collection.date == max(collection.dates),]

            first.freqs <- hist(
                first.collected$days.before.art,
                breaks=breakpoints,
                plot=FALSE
            )
            last.freqs <- hist(
                last.collected$days.before.art,
                breaks=breakpoints,
                plot=FALSE
            )

            # Sanity check: the sum of these frequencies should equal actual.freqs$count.
            if (any(actual.freqs$counts != first.freqs$counts + last.freqs$counts)) {
                cat(
                    "Warning: counts don't add up for case vl.info=",
                    vl.info,
                    ", pid=",
                    pid,
                    "\n",
                    sep=""
                )
            }

            total <- sum(first.freqs$counts, last.freqs$counts)

            # First plot the earlier ones as with the other plots.
            lines(
                seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
                first.freqs$counts / total,
                type="h",
                lwd=20,
                col=art.1.colour,
                lend="butt"
            )

            # Then plot the later ones on top.
            segments(
                seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
                y0=first.freqs$counts / total,
                y1=(first.freqs$counts + last.freqs$counts) / total,
                lwd=20,
                col=art.2.colour,
                lend="butt"
            )
        }

        # Perform the GLM fit by first re-binning all the data into *days*.
        by.day.breakpoints <- seq(0, max(curr.data$days.before.art) + 1)
        freqs.by.day <- hist(
            curr.data$days.before.art,
            breaks=by.day.breakpoints,
            plot=FALSE,
            right=FALSE
        )
        regression.frame <- data.frame(
            x=0:(length(freqs.by.day$counts) - 1),
            y=freqs.by.day$counts
        )
        decay.rate.regression <- glm(
            y ~ x,
            family=poisson,
            data=regression.frame
        )
        all.regressions[[pid]][[col.date]] <- decay.rate.regression

        # Get the fitted values at 0.1-year intervals from 0 to the number of years
        # of the infection (rounded up).
        x.fit.values <- data.frame(x=seq(0, length(actual.freqs$counts), by=0.1))
        x.fit.values.days <- x.fit.values * 365
        fit <- predict(decay.rate.regression, newdata=x.fit.values.days, se.fit=TRUE)

        predictions <- exp(fit$fit)
        upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
        lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)


        # Overlay the GLM fitted value and error bars.
        lines(
            rev(x.fit.values$x),
            365 * predictions / total.count,
            lwd=4,
            col="red"
        )
        lines(
            rev(x.fit.values$x),
            365 * upper.bounds / total.count,
            lty=2,
            lwd=2,
            col="red"
        )
        lines(
            rev(x.fit.values$x),
            365 * lower.bounds / total.count,
            lty=2,
            lwd=2,
            col="red"
        )

        # Add text listing the half life.
        decay.rate.summary <- summary(decay.rate.regression)
        x.coef <- decay.rate.summary$coefficients[2, 1]
        x.se <- decay.rate.summary$coefficients[2, 2]
        
        # t_{1/2} = - log(2) / x.coef
        half.life <- (- log(2) / x.coef) / 365
        half.life.upper <- Inf
        half.life.upper.str <- "\u221e"
        if (x.coef + 1.96 * x.se < 0) {
            half.life.upper <- (- log(2) / (x.coef + 1.96 * x.se)) / 365
            half.life.upper.str <- round(half.life.upper, digits=2)
        }
        half.life.lower <- (- log(2) / (x.coef - 1.96 * x.se)) / 365

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

        text(
            x=0,
            y=max.y * 0.96,
            label=pid,
            pos=4,
            cex=2
        )

        text(
            x=0,
            y=max.y * 0.87,
            label=substitute(
                paste(
                    t[1/2],
                    " = ",
                    hl,
                    " yr",
                    sep=""
                ),
                list(hl=round(half.life, digits=2))
            ),
            pos=4,
            cex=2
        )

        text(
            x=0,
            y=max.y * 0.795,
            label=paste(
                "(95% CI (",
                round(half.life.lower, digits=2),
                ", ",
                half.life.upper.str,  # we either already rounded it, or it's the infinity symbol
                "))",
                sep=""
            ),
            pos=4,
            cex=2
        )

        dev.off()
    }
}


####
# Make half-life plots with the half-lives computed with this method.

# Our colour scheme for reservoir sampling points.
alpha <- as.integer(0.7 * 255)
art.1.colour.alpha <- col2rgb(art.1.colour)
art.1.colour.alpha <- rgb(
    art.1.colour.alpha[1], 
    art.1.colour.alpha[2], 
    art.1.colour.alpha[3], 
    alpha, 
    maxColorValue=255
)
art.2.colour.alpha <- col2rgb(art.2.colour)
art.2.colour.alpha <- rgb(
    art.2.colour.alpha[1], 
    art.2.colour.alpha[2], 
    art.2.colour.alpha[3], 
    alpha, 
    maxColorValue=255
)


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

# Define the upper extent of the y-axis.
max.y <- 20

pdf("half_lives_no_binning.pdf")
par(mar=c(11.5, 6.5, 2, 2) + 0.1)
plot(
    x=c(1, nrow(plot.1.data) + 4),
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
    at=c(1:nrow(plot.1.data), nrow(plot.1.data) + 2, nrow(plot.1.data) + 4),
    labels=c(as.character(plot.1.data$pid), "reservoir", "proviral DNA"),
    cex.axis=2,
    las=2
)
axis(2, cex.axis=2)


# First plot the error bars (so that they'll be behind the points).  
# We do the top and bottom halves separately.
arrows(
    x0=1:nrow(plot.1.data),
    y0=plot.1.data$half.life,
    y1=plot.1.data$half.life.lower,
    angle=90,
    length=0.1,
    code=2,
    lwd=1
)
inf.upper.y <- 17  # how high the infinite upper bound arrows should be placed
upper.bound.inf <- is.infinite(plot.1.data$half.life.upper)
arrows(
    x0=which(!upper.bound.inf),
    y0=plot.1.data$half.life[!upper.bound.inf],
    y1=plot.1.data$half.life.upper[!upper.bound.inf],
    angle=90,
    length=0.1,
    code=2,
    lwd=1
)
arrows(
    x0=which(upper.bound.inf),
    y0=plot.1.data$half.life[upper.bound.inf],
    y1=inf.upper.y,
    angle=30,
    length=0.25,
    code=2,
    lwd=1
)

# We need to plot the PIDs with two sampling points separately so we can use
# a different symbol (a circle with art.1.colour as the outline and art.2.colour
# as the background).
two.sampling.indices <- plot.1.data$pid %in% two.sampling.pids
points(
    x=which(!two.sampling.indices),
    y=plot.1.data$half.life[!two.sampling.indices],
    cex=3,
    pch=19,
    col=art.1.colour.alpha
)
points(
    x=which(two.sampling.indices),
    y=plot.1.data$half.life[two.sampling.indices],
    cex=3,
    lwd=5,
    pch=21,
    col=art.1.colour.alpha,
    bg=art.2.colour.alpha
)

# Now plot the QVOA and DNA half lives.
points(
    x=c(nrow(plot.1.data) + 2, nrow(plot.1.data) + 4),
    y=c(qvoa.half.life, dna.half.life),
    cex=4,
    pch=18
)
arrows(
    x0=nrow(plot.1.data) + 2,
    y0=qvoa.half.life.lower,
    y1=qvoa.half.life.upper,
    angle=90,
    length=0.1,
    code=3,
    lwd=1
)
arrows(
    x0=nrow(plot.1.data) + 4,
    y0=dna.half.life,
    y1=dna.half.life.lower,
    angle=90,
    length=0.1,
    code=2,
    lwd=1
)
arrows(
    x0=nrow(plot.1.data) + 4,
    y0=dna.half.life,
    y1=inf.upper.y,
    angle=30,
    length=0.25,
    code=2,
    lwd=1
)

abline(v=nrow(plot.1.data) + 1, lty="dashed")

text(
    x=6,
    y=19,
    labels="untreated infection",
    cex=1.75
)

text(
    x=nrow(plot.1.data) + 2.825,
    y=19,
    labels="on ART",
    cex=1.75
)

dev.off()


# Plot 2: plot the half-lives for those individuals with two sampling
# points both separately and with their integration data combined.
plot.2.rows <- half.lives$pid %in% two.sampling.pids
plot.2.data <- half.lives[plot.2.rows,]

first.sampling <- c(1, 4, 7, 10)
second.sampling <- c(2, 5, 8, 11)
combined.sampling <- c(3, 6, 9, 12)
qvoa.x <- 17
dna.x <- 19

row.to.x <- c(1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15)

max.y <- 20

pdf("half_lives_two_sampling_points_no_binning.pdf")
par(mar=c(11.5, 6.5, 2, 2) + 0.1)
plot(
    x=c(1, dna.x),
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
    # Only label the middle of each triplet of columns.
    at=c(row.to.x[second.sampling], qvoa.x, dna.x),
    labels=c(as.character(plot.2.data$pid[second.sampling]), "reservoir", "proviral DNA"),
    cex.axis=2,
    las=2
)
axis(2, cex.axis=2)

abline(v=c(4, 8, 12), lty="dashed", col="grey")
abline(v=qvoa.x - 1, lty="dashed")

legend.loc <- legend(
    "topleft",
    legend=c("sampling 1", "sampling 2", "combined"),
    col=c(art.1.colour.alpha, art.2.colour.alpha, art.1.colour.alpha),
    lty=NA,
    lwd=c(1, 1, 5),
    pch=c(19, 19, 21),
    pt.cex=3,
    pt.bg=c(art.1.colour.alpha, art.2.colour.alpha, art.2.colour.alpha),
    cex=1.5,
    bg="white"
)

# First, plot the error bars (so the circles will be on top of them).
upper.bound.inf <- is.infinite(plot.2.data$half.life.upper)
arrows(
    x0=row.to.x[which(!upper.bound.inf)],
    y0=plot.2.data$half.life.lower[!upper.bound.inf],
    y1=plot.2.data$half.life.upper[!upper.bound.inf],
    angle=90,
    length=0.1,
    code=3,
    lwd=1
)

# For the bars that extend to infinity, we plot them with
# the upper-bound marked with an arrow.  We need special handling 
# for the arrows that are near the legend.

arrows(
    x0=row.to.x[which(upper.bound.inf)],
    y0=plot.2.data$half.life[upper.bound.inf],
    y1=plot.2.data$half.life.lower[upper.bound.inf],
    angle=90,
    length=0.1,
    code=2,
    lwd=1
)

inf.arrow.y <- 16  # how high on the plot the arrow is plotted
inf.upper.bound.x <- row.to.x[which(upper.bound.inf)]
inf.upper.bound.y <- sapply(
    inf.upper.bound.x,
    function (x) {
        if (x <= legend.loc$rect$left + legend.loc$rect$w) {
            return(legend.loc$rect$top - legend.loc$rect$h - 0.5)
        }
        return(inf.arrow.y)
    }
)

arrows(
    x0=inf.upper.bound.x,
    y0=plot.2.data$half.life[upper.bound.inf],
    y1=inf.upper.bound.y,
    angle=30,
    length=0.25,
    code=2,
    lwd=1
)

points(
    x=row.to.x[first.sampling],
    y=plot.2.data$half.life[first.sampling],
    cex=3,
    pch=19,
    col=art.1.colour.alpha
)

points(
    x=row.to.x[second.sampling],
    y=plot.2.data$half.life[second.sampling],
    cex=3,
    pch=19,
    col=art.2.colour.alpha
)

points(
    x=row.to.x[combined.sampling],
    y=plot.2.data$half.life[combined.sampling],
    cex=3,
    lwd=5,
    pch=21,
    col=art.1.colour.alpha,
    bg=art.2.colour.alpha
)

# Now plot the QVOA and DNA half lives.
points(
    x=c(qvoa.x, dna.x),
    y=c(qvoa.half.life, dna.half.life),
    cex=4,
    pch=18
)
arrows(
    x0=qvoa.x,
    y0=qvoa.half.life.lower,
    y1=qvoa.half.life.upper,
    angle=90,
    length=0.1,
    code=3,
    lwd=1
)
arrows(
    x0=dna.x,
    y0=dna.half.life,
    y1=dna.half.life.lower,
    angle=90,
    length=0.1,
    code=2,
    lwd=1
)
arrows(
    x0=dna.x,
    y0=dna.half.life,
    y1=inf.arrow.y,
    angle=30,
    length=0.25,
    code=2,
    lwd=1
)

text(
    x=12,
    y=18.5,
    labels="untreated\ninfection",
    cex=2
)

text(
    x=qvoa.x + 1,
    y=18.5,
    labels="on\nART",
    cex=2
)

dev.off()
