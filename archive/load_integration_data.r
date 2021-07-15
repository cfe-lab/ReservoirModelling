brooks.data <- read.csv("../data/BrooksData.csv")[,1:8]

names(brooks.data) <- c(
    "pid",
    "est.infection.date",
    "art.start.date",
    "seq.id",
    "collection.date",
    "integration.date.est",
    "integration.date.lower",
    "integration.date.upper"
)

for (col.idx in c(2, 3, 5, 6, 7, 8)) {
    brooks.data[[col.idx]] <- strptime(brooks.data[[col.idx]], "%d-%b-%y")
}

# Get the gap between infection date and ART initiation.
brooks.data$untreated.period <-
    as.numeric(brooks.data$art.start.date - brooks.data$est.infection.date, units="days")

# Get the number of days prior to ART initiation that each sequence is estimated
# to have been introduced to the reservoir.
brooks.data$days.before.art.raw <-
    as.numeric(brooks.data$art.start.date - brooks.data$integration.date.est, units="days")

brooks.data$days.before.art <-
    sapply(
        brooks.data$days.before.art.raw,
        function (x) {
            if (x >= 0) {
                return(x)
            }
            return(0)
        }
    )

integration.data <- brooks.data[, c(1, 5, 9, 11)]
