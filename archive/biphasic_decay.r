source("reservoir_helpers.r")

bin.size = 365
# Perform the ODE simulation over a 10 year period.
undecayed.result <- undecayed.reservoir.distribution(bin.size * 10)

# Decay the reservoir using a half life of 180 days.
decay.phase.one <- decay.distribution(undecayed.result$bin.freqs, 180, bin.size)

# Then decay each bin using a half life of 140 months for, say, 5 years.
decay.phase.two.freqs <- sapply(
    1:length(decay.phase.one$bin.freqs),
    function(idx) {
        reservoir.decay(
            decay.phase.one$bin.freqs[idx],
            140 * bin.size,  # this is the half-life
            5 * 12 * bin.size  # decay it for this long
        )
    }
)
decay.phase.two.dist <- decay.phase.two.freqs / sum(decay.phase.two.freqs)

x.labels <- length(undecayed.result$bin.freqs) - 1:length(undecayed.result$bin.freqs) + 1

pdf("monophasic_vs_biphasic_decay.pdf")
barplot(
    rbind(decay.phase.one$bin.dist, decay.phase.two.dist),
    beside=TRUE,
    names=x.labels,
    main="Reservoir composition",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir",
    legend.text=c("monophasic", "biphasic"),
    args.legend=list(x="topleft")
)
dev.off()
