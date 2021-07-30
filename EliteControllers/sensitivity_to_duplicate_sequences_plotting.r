# This is commented out as we will often have already loaded it when we run this script.
# load("sensitivity_to_duplicate_sequences.RData")
library(plotrix)

source("plot_helpers.r")

subjects <- names(all.log.likelihoods)
regime <- "Miura"

# Make plots of the decay rate estimates.
gaps <- list()
x.tick.marks <- list()
for (subject in subjects) {
    gaps[[subject]] <- list()
    x.tick.marks[[subject]] <- list()
    gaps[[subject]][[regime]] <- c(7500, 49500)
    x.tick.marks[[subject]][[regime]] <- c(seq(0, 7000, 1000), seq(50000, max(possible.half.lives), 1000))
}

# Special handling for p2-Miura.
gaps$p2$Miura <- c(6500, 43500)
x.tick.marks$p2$Miura <- c(seq(0, 6000, 1000), seq(44000, max(possible.half.lives), 1000))


for (subject in subjects) {
    mle.row <- mles[mles$subject == subject & mles$regime == regime,]
    pdf(paste("lls_", subject, "_", regime, "_with_duplicates_known_vl.pdf", sep=""))
    ll.plot(
        all.log.likelihoods[[subject]][[regime]],
        possible.half.lives,
        mle.row$mle,
        mle.row$lower.bound,
        mle.row$upper.bound,
        gaps[[subject]][[regime]],
        x.tick.marks[[subject]][[regime]]
    )
    dev.off()
}


# Plots of reservoir composition.
for (subject in subjects) {
    curr.data <- integration.data[[subject]]

    days.pre.therapy <- 
        as.numeric(
            vl.data[[subject]]$art.initiation - vl.data[[subject]]$infection.date,
            units="days"
        )

    mle.row <- mles[mles$subject == subject & mles$regime == regime,]

    # This is our criteria for whether our model is supported or whether a 
    # "no-decay" model is a better fit.
    bayes.factor <- bayes.factors$bayes.factor[
        bayes.factors$subject == subject 
        & bayes.factors$regime == regime
    ]

    x.label <- "Year prior to ART initiation"
    legend.location <- "topleft"
    if (subject == "p3") {
        x.label <- "Year prior to last viremic episode"
        legend.location <- "topright"
    }

    pdf(paste("composition_", subject, "_", regime, "_with_duplicates_known_vl.pdf", sep=""))
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
