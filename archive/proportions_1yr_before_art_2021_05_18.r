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

prop.1yr <- NULL
for (pid in unique(brooks.data$pid)) {
    curr.data <- brooks.data[brooks.data$pid == pid,]
    total.seqs <- nrow(curr.data)
    yr.before.art <- sum(curr.data$days.before.art <= 365)
    prop.1yr <- rbind(
        prop.1yr,
        data.frame(
            pid=pid,
            total=total.seqs,
            yr.before.art=yr.before.art,
            proportion=yr.before.art / total.seqs
        )
    )
}

mean.1yr <- mean(prop.1yr$proportion)
median.1yr <- median(prop.1yr$proportion)
overall.proportion <- sum(brooks.data$days.before.art <= 365) / nrow(brooks.data)

prop.1yr.summary <- data.frame(
    mean=mean.1yr,
    median=median.1yr,
    overall=overall.proportion
)

write.csv(prop.1yr, "proportions_1yr.csv")
write.csv(prop.1yr.summary, "proportions_1yr_summary.csv")
