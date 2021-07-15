model.free.regression <- function (days.before.art, breakpoints) {
    actual.freqs <- hist(
        curr.data$days.before.art,
        breaks=breakpoints,
        plot=FALSE
    )        
    regression.frame <- data.frame(
        x=1:length(actual.freqs$counts),
        y=actual.freqs$counts
    )
    return(
        list(
            actual.freqs=actual.freqs,
            fit=glm(
                y ~ x,
                family=poisson,
                data=regression.frame
            )
        )
    )
}


model.free.half.lives <- function (
    model.free.fit,
    glm.scale=1  # set this to 365 if the fit was binning by day
) {
    decay.rate.summary <- summary(model.free.fit)
    x.coef <- decay.rate.summary$coefficients[2, 1]
    x.se <- decay.rate.summary$coefficients[2, 2]
    
    # t_{1/2} = - log(2) / x.coef
    half.life <- (- log(2) / x.coef) / glm.scale
    half.life.upper <- Inf
    if (x.coef + 1.96 * x.se < 0) {
        half.life.upper <- (- log(2) / (x.coef + 1.96 * x.se)) / glm.scale
    }
    half.life.lower <- (- log(2) / (x.coef - 1.96 * x.se)) / glm.scale

    return(
        data.frame(
            half.life=half.life,
            lower.bound=half.life.lower,
            upper.bound=half.life.upper
        )
    )
}


# Our colour scheme for the "bar" parts of these graphs.
default.art.1.colour <- "#5581B0"
default.art.2.colour <- "#74DCD0"

model.free.decay.plot <- function(
    decay.rate.regression,
    breakpoints,
    integration.data,
    pid,
    glm.scale=1,  # set this to 365 if the fit was using data binned by day
    art.1.colour=default.art.1.colour,
    art.2.colour=default.art.2.colour
) {
    actual.freqs <- hist(
        integration.data$days.before.art,
        breaks=breakpoints,
        plot=FALSE
    )

    total.count <- sum(actual.freqs$counts)
    max.y <- 1

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

    # Plot the empirical distribution as a stacked bar graph.  
    # There are at most two collection dates.
    collection.dates <- unique(integration.data$collection.date)
    if (length(collection.dates) == 1) {
        lines(
            seq(length(actual.freqs$counts), 1, by=-1) - 0.5,
            actual.freqs$counts / total.count,
            type="h",
            lwd=20,
            col=art.1.colour,
            lend="butt"
        )
    } else {
        first.collected <- integration.data[integration.data$collection.date == min(collection.dates),]
        last.collected <- integration.data[integration.data$collection.date == max(collection.dates),]

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
                "Warning: counts don't add up for pid=",
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

    # Overlay the GLM fitted value and error bars.
    x.fit.values <- data.frame(x=seq(0.5, length(actual.freqs$counts) + 0.5, by=0.1))
    fit <- predict(decay.rate.regression, newdata=x.fit.values * glm.scale, se.fit=TRUE)

    predictions <- exp(fit$fit)
    upper.bounds <- exp(fit$fit + 1.96 * fit$se.fit)
    lower.bounds <- exp(fit$fit - 1.96 * fit$se.fit)

    lines(
        rev(x.fit.values$x) - 0.5,
        glm.scale * predictions / total.count,
        lwd=4,
        col="red"
    )
    lines(
        rev(x.fit.values$x) - 0.5,
        glm.scale * upper.bounds / total.count,
        lty=2,
        lwd=2,
        col="red"
    )
    lines(
        rev(x.fit.values$x) - 0.5,
        glm.scale * lower.bounds / total.count,
        lty=2,
        lwd=2,
        col="red"
    )

    # Add text listing the half life.
    hl.df <- model.free.half.lives(decay.rate.regression, glm.scale=glm.scale)

    half.life <- hl.df$half.life
    half.life.upper <- "\u221e"
    if (is.finite(hl.df$upper.bound)) {
        half.life.upper <- round(hl.df$upper.bound, digits=2)
    }
    half.life.lower <- hl.df$lower.bound

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
            half.life.upper,  # we either already rounded it, or it's the infinity symbol
            "))",
            sep=""
        ),
        pos=4,
        cex=2
    )
}
