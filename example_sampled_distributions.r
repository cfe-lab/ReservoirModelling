# Make "barplots" with none of the lines for:
# - N133M combined 
# - Z634F combined
# That is, these are like the composition plots but with only the 
# bars representing the actual empirical distribution of the reservoir,
# without any of the "modelled" stuff.

# This is commented out because we will often already have this loaded.
source("ode_based_helpers.r")
integration.data <- load.integration.data()


# Make "barplots" with none of the lines for N133M combined and Z634F combined.
art.1.colour <- "#5581B0"
art.2.colour <- "#74DCD0"
vl.info <- "typical"
for (pid in c("N133M", "Z634F")) {
    curr.data <- integration.data[integration.data$pid == pid,]
    days.pre.therapy <- curr.data$untreated.period[1]

    breakpoints = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=365)
    if (!(max(curr.data$days.before.art) %in% breakpoints)) {
        breakpoints = c(breakpoints, breakpoints[length(breakpoints)] + 365)
    }
    actual.freqs <- hist(
        curr.data$days.before.art,
        breaks=breakpoints,
        plot=FALSE
    )
    emp.dist <- actual.freqs$counts / sum(actual.freqs$counts)

    max.y <- 1

    pdf(paste("sampled_", pid, ".pdf", sep=""))
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

    # Plot the empirical distribution as a stacked bar graph.  
    # These samples both have  two collection dates.
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

    legend.labels <- c(
        "sampling 1",
        "sampling 2"
    )
    legend.colours <- c(art.1.colour, art.2.colour)
    legend.ltys <- c("solid", "solid")
    legend.lwds <- c(20, 20)

    legend.location <- legend(
        "topleft",
        legend=legend.labels,
        col=legend.colours,
        lty=legend.ltys,
        lwd=legend.lwds,
        bty="n",
        cex=1.5
    )

    text(
        x=length(emp.dist) * 0.8,
        y=0.95,
        labels=pid,
        cex=2
    )
    dev.off()
}
