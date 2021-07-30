# This is commented out as we will often have already loaded it when we run this script.
# load("p3_no_special_handling_ode.RData")
# source("p3_no_special_handling_include_duplicates_processing.r")
library(plotrix)

source("plot_helpers.r")


# Make plots of the decay rate estimates.
gaps <- list()
x.tick.marks <- list()
for (regime in regimes) {
    gaps[[regime]] <- c(7500, 49500)
    x.tick.marks[[regime]] <- c(seq(0, 7000, 1000), seq(50000, max(possible.half.lives), 1000))
}


subject <- "p3"
for (regime in regimes) {
    mle.row <- mles[mles$subject == subject & mles$regime == regime,]

    pdf(paste("lls_", subject, "_", regime, "_known_vl_no_special_handling_include_duplicates.pdf", sep=""))
    ll.plot(
        all.log.likelihoods[[subject]][[regime]],
        possible.half.lives,
        mle.row$mle,
        mle.row$lower.bound,
        mle.row$upper.bound,
        gaps[[regime]],
        x.tick.marks[[regime]]
    )
    dev.off()
}


# Plots of reservoir composition.
curr.data <- integration.data[[subject]]

days.pre.therapy <- 
    as.numeric(
        vl.data[[subject]]$art.initiation - vl.data[[subject]]$infection.date,
        units="days"
    )

for (regime in regimes) {
    mle.row <- mles[mles$subject == subject & mles$regime == regime,]

    # This is our criteria for whether our model is supported or whether a 
    # "no-decay" model is a better fit.
    bayes.factor <- bayes.factors$bayes.factor[
        bayes.factors$subject == subject 
        & bayes.factors$regime == regime
    ]

    # Recall that there's no special handling for P3.
    x.label <- "Year prior to ART initiation"
    legend.location <- "topleft"

    pdf(paste("composition_", subject, "_", regime, "_known_vl_no_special_handling_include_duplicates.pdf", sep=""))
    composition.plot(
        ode.solutions$bin.365[[subject]][[regime]]$bin.freqs,
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
