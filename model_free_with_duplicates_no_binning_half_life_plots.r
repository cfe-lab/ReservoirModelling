# This is akin to model_free_no_binning_half_life_plots.r but with
# the integration data that includes duplicate proviruses.

source("ode_based_helpers.r")
integration.data <- load.combined.integration.data(
    load.integration.data(),
    load.integration.data.with.duplicates()
)

source("model_free_helpers.r")
source("half_life_plot_helpers.r")

# This will be doubly-indexed by:
#  - PID
#  - collection date (could also be "combined")
all.regressions <- list()
half.lives <- NULL


for (pid in unique(integration.data$pid)) {
    all.regressions[[pid]] <- list()

    curr.pid.data <- integration.data[integration.data$pid == pid,]
    days.pre.therapy <- curr.pid.data$untreated.period[1]

    col.dates <- as.character(unique(curr.pid.data$collection.date))
    if (length(col.dates) > 1) {
        col.dates <- c(col.dates, "combined")
    }

    for (col.date in col.dates) {
        curr.data <- curr.pid.data
        if (col.date != "combined") {
            curr.data <- curr.pid.data[as.character(curr.pid.data$collection.date) == col.date,]
        }

        # The GLM fit will use the data binned by day.
        bp.days <- seq(0, max(curr.data$days.before.art) + 1)

        curr.mfr <- model.free.regression(
            curr.data$days.before.art,
            bp.days
        )
        all.regressions[[pid]][[col.date]] <- curr.mfr

        curr.hl <- model.free.half.lives(curr.mfr$fit, glm.scale=365)
        half.lives <- rbind(
            half.lives,
            cbind(
                data.frame(
                    pid=pid,
                    col.date=col.date
                ),
                curr.hl
            )
        )

        # We still need to bin the data by year for the plot, so we
        # compute those breakpoints here.
        bp.years = seq(0, max(c(curr.data$days.before.art, days.pre.therapy)), by=365)
        if (!(max(curr.data$days.before.art) %in% bp.years)) {
            bp.years = c(bp.years, bp.years[length(bp.years)] + 365)
        }

        cairo_pdf(
            paste0(
                "model_free_with_duplicates_no_binning_decay_rate_", 
                pid,
                "_",
                col.date,
                ".pdf"
            )
        )
        model.free.decay.plot(
            curr.mfr$fit,
            bp.years,
            curr.data,
            pid,
            glm.scale=365
        )
        dev.off()
    }
}


# Plot the half-lives for all individuals, only including the "combined"
# sampling points for those with multiple sampling timepoints.
pdf("model_free_with_duplicates_no_binning_half_lives.pdf")
plot.all.individuals.half.lives(half.lives)
dev.off()


# Plot the half-lives for those individuals with two sampling
# points both separately and with their integration data combined.
pdf("model_free_with_duplicates_no_binning_half_lives_two_sampling_points.pdf")
plot.multiple.timepoints.half.lives(half.lives)
dev.off()
