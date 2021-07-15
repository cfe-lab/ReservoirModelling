# Our colour scheme for the "bar" parts of these graphs.
default.art.1.colour <- "#5581B0"
default.art.2.colour <- "#74DCD0"
default.alpha <- as.integer(0.7 * 255)


# A function to build all of our half-life plots (the scatterplots with error bars).
plot.half.lives <- function(
    half.lives,  # columns: label, col.date ("first", "second", or "combined"), half.life, upper.bound, lower.bound
    separate.after.row,  # a vector of row indices to place a vertical separator after
    show.legend=FALSE,
    art.1.colour=default.art.1.colour,
    art.2.colour=default.art.2.colour,
    alpha=default.alpha
) {
    # Our colour scheme for reservoir sampling points.
    alpha <- as.integer(0.7 * 255)
    art.1.colour <- col2rgb(art.1.colour)
    art.1.colour <- rgb(art.1.colour[1], art.1.colour[2], art.1.colour[3], alpha, maxColorValue=255)
    art.2.colour <- col2rgb(art.2.colour)
    art.2.colour <- rgb(art.2.colour[1], art.2.colour[2], art.2.colour[3], alpha, maxColorValue=255)

    # We'll plot all the half-lives alongside the QVOA and DNA half-lives.
    # The numerators below are all in months so we divide by 12 to get years.
    qvoa.half.life <- 44.2 / 12
    qvoa.half.life.upper <- 114.5 / 12
    qvoa.half.life.lower <- 27.4 / 12
    dna.half.life <- 140.4 / 12
    dna.half.life.upper <- Inf
    dna.half.life.lower <- 75.6 / 12

    # We map our row indices in half.lives to x-axis coordinates, factoring
    # in where the vertical partitions go.
    row.to.x <- NULL
    vertical.partitions <- NULL
    curr.row.idx <- 1
    for (row.idx in 1:nrow(half.lives)) {
        row.to.x <- c(row.to.x, curr.row.idx)
        curr.row.idx <- curr.row.idx + 1
        if (row.idx %in% separate.after.row) {
            vertical.partitions <- c(vertical.partitions, curr.row.idx)
            curr.row.idx <- curr.row.idx + 1
        }
    }

    # Define the upper extent of the y-axis, and the y-value to which "infinite" arrows
    # will extend.
    max.y <- 20
    inf.arrow.y <- 17

    par(mar=c(11.5, 6.5, 2, 2) + 0.1)
    plot(
        x=c(1, max(row.to.x) + 4),  # includes space for the QVOA and DNA points
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

    rows.with.labels <- half.lives$label != ""
    x.labels <- half.lives$label[rows.with.labels]
    x.label.positions <- row.to.x[rows.with.labels]

    axis(
        1,  # this is the x axis
        at=c(x.label.positions, max(row.to.x) + 2, max(row.to.x) + 4),
        labels=c(x.labels, "reservoir", "proviral DNA"),
        cex.axis=2,
        las=2
    )
    axis(2, cex.axis=2)

    # First, plot the vertical partitions, so that the legend can go on top of them.
    for (x.partition in vertical.partitions) {
        abline(v=x.partition, lty="dashed", col="grey")
    }
    abline(v=max(row.to.x) + 1, lty="dashed")

    left.top.untreated.area <- 1
    right.top.untreated.area <- max(row.to.x) + 1
    # Next, plot the legend, so we can store its location.
    if (show.legend) {
        legend.loc <- legend(
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
        left.top.untreated.area <- legend.loc$rect$left + legend.loc$rect$w
    }

    # First plot the error bars (so that they'll be behind the points).  
    # We do the top and bottom halves separately.
    arrows(
        x0=row.to.x[1:nrow(half.lives)],
        y0=half.lives$half.life,
        y1=half.lives$lower.bound,
        angle=90,
        length=0.1,
        code=2,
        lwd=1
    )

    # There are three cases for the upper bounds:
    # - finite and easily plottable
    # - infinite
    # - finite but too big to plot
    for (idx in 1:nrow(half.lives)) {
        x.coord <- row.to.x[idx]
        raw.upper.bound <- half.lives$upper.bound[idx]
        hl <- half.lives$half.life[idx]
        y.extent <- inf.arrow.y
        if (show.legend && x.coord <= legend.loc$rect$left + legend.loc$rect$w) {
            y.extent <- legend.loc$rect$top - legend.loc$rect$h - 0.5
        }

        if (is.finite(raw.upper.bound)) {
            if (raw.upper.bound <= y.extent) {
                # This is the most normal case.
                arrows(
                    x0=x.coord,
                    y0=hl,
                    y1=raw.upper.bound,
                    angle=90,
                    length=0.1,
                    code=2,
                    lwd=1
                )
            } else {
                # Instead of an arrow, we draw a solid line to y.extent - 3,
                # a dotted line to y.extent - 2, and then put the actual upper bound in writing.
                segments(
                    x0=x.coord,
                    y0=hl,
                    y1=y.extent - 3
                )
                segments(
                    x0=x.coord,
                    y0=y.extent - 3,
                    y1=y.extent - 2,
                    lty="dotted"
                )
                text(
                    x.coord,
                    y.extent - 1.25,
                    labels=round(raw.upper.bound, digits=2),
                    cex=1.25
                )
            }
        } else {
            # Draw a "pointing" arrow that extends up to y.extent.
            arrows(
                x0=x.coord,
                y0=hl,
                y1=y.extent,
                angle=30,
                length=0.125,
                code=2,
                lwd=1
            )
        }
    }

    # Plot points for each row, colouring them by their collection date.
    first.sampling <- half.lives$col.date == "first"
    second.sampling <- half.lives$col.date == "second"
    combined.sampling <- half.lives$col.date == "combined"
    points(
        x=row.to.x[first.sampling],
        y=half.lives$half.life[first.sampling],
        cex=3,
        pch=19,
        col=art.1.colour
    )
    points(
        x=row.to.x[second.sampling],
        y=half.lives$half.life[second.sampling],
        cex=3,
        pch=19,
        col=art.2.colour
    )
    points(
        x=row.to.x[combined.sampling],
        y=half.lives$half.life[combined.sampling],
        cex=3,
        lwd=5,
        pch=21,
        col=art.1.colour,
        bg=art.2.colour
    )

    # Now plot the QVOA and DNA half lives.
    qvoa.x <- max(row.to.x) + 2
    dna.x <- max(row.to.x) + 4
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
        length=0.125,
        code=2,
        lwd=1
    )

    text(
        x=(left.top.untreated.area + right.top.untreated.area) / 2,
        y=18.5,
        # labels="untreated infection",
        labels="untreated",
        cex=2
    )

    text(
        x=max(row.to.x) + 3,
        y=18.5,
        # labels="on ART",
        labels="ART",
        cex=2
    )
}


identify.duplicate.pids <- function(pids) {
    duplicate.pids <- NULL
    for (pid in unique(pids)) {
        if (sum(pids == pid) > 1) {
            duplicate.pids <- c(duplicate.pids, pid)
        }
    }
    return(duplicate.pids)
}


# A helper that produces the plots of half-lives for everyone, where
# those individuals with more than one sampling timepoint only have
# their "combined" entries plotted.
plot.all.individuals.half.lives <- function(half.lives) {
    two.sampling.pids <- identify.duplicate.pids(half.lives$pid)

    all.individuals.rows <- !(half.lives$pid %in% two.sampling.pids) | (half.lives$col.date == "combined")
    all.individuals.data <- half.lives[all.individuals.rows,]
    all.individuals.data$col.date <- sapply(
        all.individuals.data$col.date,
        function (x) {
            if (x == "combined") {
                return("combined")
            }
            return("first")
        }
    )
    all.individuals.data$label <- all.individuals.data$pid
    plot.half.lives(all.individuals.data, NULL)
}


# A helper to produce the plots of half-lives for those individuals
# with two sampling timepoints.  We assume that the entries in half.lives
# are sorted by PID and then by collection date.
plot.multiple.timepoints.half.lives <- function(half.lives) {
    two.sampling.pids <- identify.duplicate.pids(half.lives$pid)  # this preserves order

    two.samplings.rows <- half.lives$pid %in% two.sampling.pids
    two.samplings.data <- half.lives[two.samplings.rows,]

    # Identify which rows are the first samplings, which rows are the second samplings,
    # and which are the combined samplings.
    first.sampling <- NULL
    second.sampling <- NULL
    combined.sampling <- NULL
    for (pid in two.sampling.pids) {
        curr.col.dates <- two.samplings.data$col.date[two.samplings.data$pid == pid]
        first.date <- min(curr.col.dates[curr.col.dates != "combined"])
        second.date <- max(curr.col.dates[curr.col.dates != "combined"])

        curr.pid.rows <- two.samplings.data$pid == pid

        first.sampling <- c(first.sampling, which(curr.pid.rows & two.samplings.data$col.date == first.date))
        second.sampling <- c(second.sampling, which(curr.pid.rows & two.samplings.data$col.date == second.date))
        combined.sampling <- c(combined.sampling, which(curr.pid.rows & two.samplings.data$col.date == "combined"))
    }

    two.samplings.data$col.date[first.sampling] <- "first"
    two.samplings.data$col.date[second.sampling] <- "second"
    # We only label the middle columns for each triplet (each representing
    # one person).
    two.samplings.data$label <- ""
    two.samplings.data$label[second.sampling] <- two.samplings.data$pid[second.sampling]

    partitions.between.individuals <- 1:(length(two.samplings.data) - 1) * 3
    plot.half.lives(
        two.samplings.data,
        partitions.between.individuals,
        show.legend=TRUE
    )
}
