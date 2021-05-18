# Make plots similar to those from the 2021-04-12 script showing the
# "typical" composition with decay, but show the compositions as
# actual separate bar plots rather than as "bar-like lines" on the 
# same plot.

# load("typical_pvl_demonstration_2021_04_07.RData")

source("reservoir_helpers.r")

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

    max.y <- 1

    pdf(paste0("composition_", num.untreated.years, "yr_untreated.pdf"), height=9.5)
    # par(mar=c(10.5, 10.5, 2, 2) + 0.1, mfrow=c(4, 1))
    par(mfrow=c(4, 1), oma=c(3, 4.5, 0, 0), xpd=NA)

    values.to.plot <- list()
    values.to.plot[[1]] <- undecayed
    values.to.plot[[2]] <- dist.140mo.decay$bin.dist
    values.to.plot[[3]] <- dist.44mo.decay$bin.dist
    values.to.plot[[4]] <- dist.180day.decay$bin.dist
    plot.labels <- c(
        "created by viral seeding: no decay",
        "with 140 mo decay",
        "with 44 mo decay",
        "with 6 mo decay"
    )
    for (idx in 1:4) {
        bar.x.values <- barplot(
            values.to.plot[[idx]],
            names.arg=seq(length(bin.freqs), 1, by=-1),
            ylim=c(0, max.y),
            xlab=NA,
            ylab=NA,
            cex.names=2,
            cex.axis=2
        )
        text(
            x=bar.x.values[trunc((num.untreated.years + 1) / 2)],
            y=0.9,
            labels=plot.labels[idx],
            cex=2
        )   
    }

    mtext(
        text="Year prior to ART initiation",
        side=1,
        outer=TRUE,
        cex=2.25
    )

    mtext(
        text="Proportion of proviruses",
        side=2,
        outer=TRUE,
        cex=3
    )

    dev.off()
}
