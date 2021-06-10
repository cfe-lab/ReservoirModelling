# For p2, we compute a likelihood ratio/Bayes factor to compare the models with
# - MLE decay
# - no decay at all

# load("full_analysis.RData")

subject <- "p2"
bayes.factors.p2 <- NULL
for (regime in regimes) {
    curr.data <- all.subjects[[subject]]$integration
    reservoir.dist <- ode.solutions$bin.30[[subject]][[regime]]

    no.decay.ll.no.factorial <- log.likelihood.no.factorial(
        curr.data$days.before.art,
        reservoir.dist$bin.dist.no.decay,
        30
    )
    no.decay.ll <- no.decay.ll.no.factorial + sum(log(1:nrow(curr.data)))

    lls <- all.log.likelihoods[[subject]][[regime]]
    max.idx <- which.max(lls)
    max.ll <- lls[max.idx]

    bayes.factor <- exp(no.decay.ll - max.ll)
    bayes.factors.p2 <- rbind(
        bayes.factors.p2,
        data.frame(
            regime=regime,
            bayes.factor=bayes.factor
        )
    )
}
