# Use a generalized linear model to find a decay rate directly from the observed
# data, without using our model.
# The details of this analysis:
#  - "standard" integration data
#  - special handling for P3's VL data (i.e. model starts from the 
#    blip rather than ART initiation) and integration dates
#  - duplicate proviruses removed

source("analysis_helpers.r")
source("plot_helpers.r")

# Read in the viral load data and proviral integration data.
read.vl <- prepare.vl.data(subjects, p3.special.handling=TRUE)
vl.data <- read.vl$vl.data
integration.data <- prepare.integration.data(
    subjects,
    vl.data,
    remove.duplicates=TRUE,
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

    # As per this link:
    # http://biometry.github.io/APES/Stats/stats23-GeneralizedLinearModels-GLM.html
    # there is *some* evidence of overdispersion in this fit, but not as decisive
    # as in the example.

    x.label="Year prior to ART initiation"
    if (subject == "p3") {
        x.label="Year prior to last viremic episode"
    }
    cairo_pdf(paste("model_free_decay_rate_", subject, ".pdf", sep=""))
    model.free.plot(
        decay.rate.regression,
        actual.freqs$counts,
        x.label=x.label
    )
    dev.off()
}
