# This is commented out as we will often have already loaded it when we run this script.
# load("p3_no_special_handling_analysis.RData")
library(plotrix)

source("plot_helpers.r")

regimes <- names(all.log.likelihoods$p3)
for (regime in regimes) {
    mle.row <- mles[mles$subject == "p3" & mles$regime == regime,]

    pdf(paste("p3_", regime, "_known_vl_no_special_handling.pdf", sep=""))
    ll.plot(
        all.log.likelihoods[["p3"]][[regime]],
        possible.half.lives,
        mle.row$mle,
        mle.row$lower.bound,
        mle.row$upper.bound,
        c(7500, 49500),
        c(seq(0, 7000, 1000), seq(50000, max(possible.half.lives), 1000))
    )
    dev.off()
}


curr.data <- all.subjects[["p3"]]$integration
days.pre.therapy <- 
    as.numeric(
        all.subjects[["p3"]]$art.initiation - all.subjects[["p3"]]$infection.date,
        units="days"
    )


# Plots of reservoir composition.
for (regime in regimes) {
    mle.row <- mles[mles$subject == "p3" & mles$regime == regime,]

    # This is our criteria for whether our model is supported or whether a 
    # "no-decay" model is a better fit.
    bayes.factor <- bayes.factors$bayes.factor[
        bayes.factors$subject == "p3" 
        & bayes.factors$regime == regime
    ]

    x.label <- "Year prior to ART initiation"  # remember: no special handling
    legend.location <- "topright"

    pdf(paste("composition_p3_", regime, "_known_vl_no_special_handling.pdf", sep=""))
    composition.plot(
        ode.solutions$bin.365[["p3"]][[regime]]$bin.freqs,
        curr.data$days.before.art,
        days.pre.therapy,
        mle.row$mle,
        mle.row$lower.bound,
        mle.row$upper.bound,
        bayes.factor=bayes.factor,
        best.fit.no.decay=bayes.factor < 1,
        x.label=x.label,
        legend.location=legend.location
    )
    dev.off()
}


# Now do the model-free decay rate calculations and plots.
breakpoints <- compute.bin.breakpoints(
    curr.data$days.before.art,
    days.pre.therapy,
    bin.size=365
)
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

cairo_pdf("model_free_decay_rate_p3_no_special_handling.pdf")
model.free.plot(
    decay.rate.regression,
    actual.freqs$counts
)
dev.off()
