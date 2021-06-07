# Jan 21, 2021: this builds on the same analysis as the 2020_12_22 script
# but saves a bit more information in an RData file for subsequent plotting.

load("elite_controllers_known_vl_2020_12_22.RData")

# Find the decay rate that maximizes the likelihood for each individual.
# NEW FOR THIS ANALYSIS:
# look up to 1800 months, which is roughly 150 years.
library(parallel)
bin.size <- 30
possible.half.lives <- (1:1800) * bin.size  # months, roughly
all.log.likelihoods <- list()  # this will be triply-indexed by pid, regime, and collection date

for (subject in subjects) {
    all.log.likelihoods[[subject]] <- list()
    curr.data <- all.subjects[[subject]]$integration

    for (regime in regimes) {
        all.log.likelihoods[[subject]][[regime]] <- list()

        reservoir.dist <- ode.solutions.bin.30[[subject]][[regime]]
        cat(
            "Calculating likelihoods on decayed distributions for subject ",
            subject,
            " (",
            regime,
            " case)....\n",
            sep=""
        )
        run.time <- system.time(
            log.likelihoods.no.factorial <- mclapply(
                possible.half.lives,
                function (x) {
                    decayed <- decay.distribution(reservoir.dist$bin.freqs, x, bin.size)
                    return(
                        log.likelihood.no.factorial(
                            curr.data$days.before.art, 
                            decayed$bin.dist, 
                            bin.size
                        )
                    )
                },
                mc.cores=3
            )
        )
        # Compute the log factorial term for the likelihood.
        all.log.likelihoods[[subject]][[regime]] <- (
            unlist(log.likelihoods.no.factorial) + 
            sum(log(1:nrow(curr.data)))
        )

        cat(
            "Likelihood calculations for subject ", 
            subject, 
            " (",
            regime,
            " case) complete.\n",
            sep=""
        )
        cat("Time elapsed:\n")
        print(run.time)
    }
}

save.image("elite_controllers_known_vl_2021_01_21.RData")
