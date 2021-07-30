# Produce plots showing the half-lives computed five ways for each participant:
# - standard integration dates, duplicate proviruses removed
# - "alternative tree" integration dates
# - nearest-neighbour-based integration dates
# - LSD integration dates
# - standard integration dates, duplicate proviruses included

source("analysis_helpers.r")
source("../half_life_plot_helpers.r")  # this is from the other project!

read.vl <- prepare.vl.data(subjects, p3.special.handling=TRUE)
vl.data <- read.vl$vl.data
p3.boundaries <- c(read.vl$p3.art.initiation, read.vl$p3.blip)

all.integration.data <- list()
all.integration.data[["standard (no duplicates)"]] <- prepare.integration.data(
    subjects,
    vl.data,
    remove.duplicates=TRUE,
    p3.boundaries=p3.boundaries
)
all.integration.data[["alternative trees"]] <- prepare.alternative.trees.integration.data(
    subjects,
    vl.data,
    p3.boundaries
)
all.integration.data[["NN"]] <- prepare.nearest.neighbour.integration.data(
    subjects,
    vl.data
)
all.integration.data[["LSD"]] <- prepare.lsd.integration.data(
    subjects,
    vl.data,
    p3.boundaries
)
all.integration.data[["standard (with duplicates)"]] <- prepare.integration.data(
    subjects,
    vl.data,
    remove.duplicates=FALSE,
    p3.boundaries=p3.boundaries
)

for (subject in subjects) {
    hl.plot.df <- NULL
    for (data.type in names(all.integration.data)) {
        curr.data <- all.integration.data[[data.type]][[subject]]
        result <- compute.model.free.estimate(
            integration.data[[subject]]$days.before.art,
            vl.data[[subject]]$infection.date,
            vl.data[[subject]]$art.initiation
        )

        hl.list <- model.free.half.life(result$decay.rate.regression)
        hl.plot.df <- rbind(
            hl.plot.df,
            data.frame(
                label=data.type,
                col.date="first",
                half.life=hl.list$half.life,
                upper.bound=hl.list$upper.bound,
                lower.bound=hl.list$lower.bound
            )
        )
    }

    # Now make the plot for this participant.
    pdf(paste0("model_free_half_lives_", subject, ".pdf"))
    plot.half.lives(
        hl.plot.df,
        separate.after.row=NULL
    )
    dev.off()
}
