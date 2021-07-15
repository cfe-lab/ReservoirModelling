# Produce plots of the likelihood function, with the MLE and its error
# bounds marked where possible.

# This is commented out because we will often already have this loaded.
# load("brooks_data_ode.RData")
source("ode_based_helpers.r")

integration.data <- load.integration.data()

# These are default values:
bin.size <- 30
possible.half.lives <- (1:200) * bin.size

# Find the decay rate that maximizes the likelihood for each individual.
estimation.results <- ode.based.likelihood.estimates(
    integration.data,
    ode.results$bin.30,
    bin.size=bin.size,
    possible.half.lives=possible.half.lives,
    parallel.cores=3
)
all.log.likelihoods <- estimation.results$all.log.likelihoods
mles <- estimation.results$mles


# Make plots of the decay rate estimates.
for (vl.info in c("known", "typical")) {
    for (pid in names(all.log.likelihoods[[vl.info]])) {
        for (col.date in names(all.log.likelihoods[[vl.info]][[pid]])) {
            lls <- all.log.likelihoods[[vl.info]][[pid]][[col.date]]

            col.date.str <- "all collection dates"
            if (col.date != "combined") {
                col.date.str <- paste("collected", col.date)
            }

            mle.row <- mles[
                mles$vl.info == vl.info &
                mles$pid == pid &
                mles$col.date == col.date,
            ]
            mle <- mle.row$mle
            lower.bound <- mle.row$lower.bound
            upper.bound <- mle.row$upper.bound

            plot.title <- paste(
                "LLs (using ",
                vl.info,
                " VL) by decay rate: ", 
                pid, 
                ", ",
                col.date.str,
                sep=""
            )

            pdf(paste(pid, "_", col.date,  "_", vl.info, "_vl.pdf", sep=""))
            plot.log.likelihoods(
                lls,
                possible.half.lives,
                plot.title,
                mle.row$mle,
                mle.row$lower.bound,
                mle.row$upper.bound
            )
            dev.off()
        }
    }
}
