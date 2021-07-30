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
all.integration.data[["standard"]] <- prepare.integration.data(
    subjects,
    vl.data,
    remove.duplicates=TRUE,
    p3.boundaries=p3.boundaries
)
all.integration.data[["alternative"]] <- prepare.alternative.trees.integration.data(
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
all.integration.data[["duplicates"]] <- prepare.integration.data(
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
            curr.data$days.before.art,
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
                lower.bound=hl.list$lower.bound,
                stringsAsFactors=FALSE
            )
        )
    }

    # Now make the plot for this participant.
    max.y <- 20
    inf.arrow.y <- 17
    text.label.y <- 18.75
    if (subject == "p2") {
        magnification <- 3
        max.y <- max.y * magnification
        inf.arrow.y <- 17 * magnification
        text.label.y <- 18.75 * magnification
    }

    pdf(paste0("model_free_half_lives_", subject, ".pdf"))
    plot.half.lives(
        hl.plot.df,
        separate.after.row=NULL,
        max.y=max.y,
        inf.arrow.y=inf.arrow.y,
        text.label.y=text.label.y
    )
    dev.off()
}
