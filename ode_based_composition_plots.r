# Plots of reservoir composition.

# This is commented out because we will often already have this loaded.
# load("brooks_data_ode.RData")
source("ode_based_helpers.r")

integration.data <- load.integration.data()
estimation.results <- ode.based.likelihood.estimates(
    integration.data,
    ode.results$bin.30
)
all.log.likelihoods <- estimation.results$all.log.likelihoods
mles <- estimation.results$mles

for (vl.info in c("known", "typical")) {
    for (pid in names(all.log.likelihoods[[vl.info]])) {
        curr.pid.data <- integration.data[integration.data$pid == pid,]
        days.pre.therapy <- curr.pid.data$untreated.period[1]

        # Remember that these are strings, and may also be the string "combined",
        # meaning that information from all collection dates was combined.
        all.collection.dates <- names(all.log.likelihoods[[vl.info]][[pid]])
        for (col.date in all.collection.dates) {
            curr.data <- curr.pid.data
            if (col.date != "combined") {
                curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
            }

            reservoir.dist <- ode.results$bin.365[[vl.info]][[pid]]

            mle.row <- mles[
                mles$vl.info == vl.info &
                mles$pid == pid &
                mles$col.date == col.date,
            ]

            pdf(paste("composition_", pid, "_", col.date,  "_", vl.info, "_vl.pdf", sep=""))
            plot.ode.based.reservoir.composition(
                curr.data$days.before.art,
                curr.data$collection.date,
                reservoir.dist$bin.freqs,
                days.pre.therapy,
                mle.row$mle,
                mle.row$lower.bound,
                mle.row$upper.bound,
                show.legend=(pid == "N133M")  # special handling for the N133M plots
            )
            dev.off()
        }
    }
}
