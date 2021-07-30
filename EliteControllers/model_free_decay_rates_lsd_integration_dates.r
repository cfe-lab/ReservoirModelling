# This is akin to model_free_decay_rates.r, but we use the integration dates
# drawn using the least-squares-dating method.
# The details of this analysis:
#  - LSD integration data
#  - special handling for P3's VL data (i.e. model starts from the 
#    blip rather than ART initiation) and integration dates

source("plot_helpers.r")
source("analysis_helpers.r")

# Read in the viral load data.
read.vl <- prepare.vl.data(subjects, p3.special.handling=TRUE)
vl.data <- read.vl$vl.data
integration.data <- prepare.lsd.integration.data(
    subjects,
    vl.data,
    p3.boundaries=c(read.vl$p3.art.initiation, read.vl$p3.blip)
)

for (subject in subjects) {
    result <- compute.model.free.estimate(
        integration.data[[subject]]$days.before.art,
        vl.data[[subject]]$infection.date,
        vl.data[[subject]]$art.initiation
    )
    actual.freqs <- result$actual.freqs
    decay.rate.regression <- result$decay.rate.regression

    x.label="Year prior to ART initiation"
    if (subject == "p3") {
        x.label="Year prior to last viremic episode"
    }
    cairo_pdf(paste("model_free_decay_rate_lsd_integration_dates_", subject, ".pdf", sep=""))
    model.free.plot(
        decay.rate.regression,
        actual.freqs$counts,
        x.label=x.label
    )
    dev.off()
}
