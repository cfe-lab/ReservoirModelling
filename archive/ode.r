# First, install deSolve.

library(deSolve)

parameters <- c(
    alpha.S=70,
    delta.S=0.2,
    beta=0.0001,
    tau=0.05,
    lambda=0.0001,
    delta.I = 0.8,
    pi=50000,
    gamma=23,
    alpha.E=0.0001,
    kappa=0.3,
    omega=1.6,
    delta.E=0.002,
    E.50=250
)

I.0 = parameters["gamma"] / parameters["pi"] * 0.03

state <- c(
    S=70 / 0.2,
    A.P=I.0 * (1 - parameters["tau"]) * (1 - parameters["lambda"]),
    A.U=I.0 * parameters["tau"] * (1 - parameters["lambda"]),
    L.P=I.0 * (1 - parameters["tau"]) * parameters["lambda"],
    L.U=I.0 * parameters["tau"] * parameters["lambda"],
    E=0.0001 / 0.002,
    # E=1,
    V=30 / 1000
)
names(state) <- c("S", "A.P", "A.U", "L.P", "L.U", "E", "V")

derivatives <- function(t, state, parameters) {
    with(
        as.list(c(state, parameters)), 
        {
            dS <- alpha.S - delta.S * S - beta * S * V
            dA.P <- (1 - lambda) * tau * beta * S * V - delta.I * A.P - kappa * A.P * E
            dA.U <- (1 - lambda) * (1 - tau) * beta * S * V - delta.I * A.U - kappa * A.U * E
            # dA.P <- tau * beta * S * V - delta.I * A.P - kappa * A.P * E
            # dA.U <- (1 - tau) * beta * S * V - delta.I * A.U - kappa * A.U * E
            dL.P <- lambda * tau * beta * S * V
            dL.U <- lambda * (1 - tau) * beta * S * V
            dE <- alpha.E + omega * (A.P + A.U) * E / (E + E.50) - delta.E * E
            dV <- pi * A.P - gamma * V - beta * S * V

            return(list(c(dS, dA.P, dA.U, dL.P, dL.U, dE, dV)))
        }
    )
}

times <- seq(0, 365 * 10, by=0.001)

out <- ode(
    y=state,
    times=times,
    func=derivatives,
    parms=parameters,
    maxsteps=25000
)

plot(out, xlab="time")

out.df <- as.data.frame(out)
pdf("modelled_vl_progression.pdf")
plot(
    out.df$time[1:180000], 
    1000 * out.df$V[1:180000],
    xlab="days post-infection",
    ylab="Viral load",
    main="Modelled progression of viral load",
    type="l",
    log="y"
)
dev.off()

year.indices <- seq(1, 365 * 10 * 1000 + 1, by=365 * 1000)

latent.by.year <- (
    out.df$L.P[year.indices[2:length(year.indices)]] - 
    out.df$L.P[year.indices[1:length(year.indices)-1]]
)

reservoir.decay <- function(size, half.life, decay.time) {
    # Parameters:
    # size is the number of particles in question.
    # half.life is the half-life in days.
    # decay.time is the amount of time the particles have aged.
    return(size * 2^(-(decay.time / half.life)))
}

x.labels <- length(latent.by.year) - 1:length(latent.by.year) + 1

pdf("no_decay.pdf")
barplot(
    latent.by.year / sum(latent.by.year),
    names.arg=x.labels,
    main="Composition of reservoir by age (no decay)",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir"
)
dev.off()

decay.44mo.hl <- sapply(
    1:length(latent.by.year),
    function (idx) {
        reservoir.decay(
            latent.by.year[idx],
            44 * 30,
            365 * (10 - idx)
        )
    }
)

pdf("44mo_decay.pdf")
barplot(
    decay.44mo.hl / sum(decay.44mo.hl),
    names.arg=x.labels,
    main="Composition of reservoir by age (44mo half-life)",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir"
)
dev.off()

decay.140mo.hl <- sapply(
    1:length(latent.by.year),
    function (idx) {
        reservoir.decay(
            latent.by.year[idx],
            140 * 30,
            365 * (10 - idx)
        )
    }
)

pdf("140mo_decay.pdf")
barplot(
    decay.140mo.hl / sum(decay.140mo.hl),
    names.arg=x.labels,
    main="Composition of reservoir by age (140mo half-life)",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir"
)
dev.off()

decay.2wk.hl <- sapply(
    1:length(latent.by.year),
    function (idx) {
        reservoir.decay(
            latent.by.year[idx],
            14,
            365 * (10 - idx)
        )
    }
)

pdf("2wk_decay.pdf")
barplot(
    decay.2wk.hl / sum(decay.2wk.hl),
    names.arg=x.labels,
    main="Composition of reservoir by age (2wk half-life)",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir"
)
dev.off()

decay.6mo.hl <- sapply(
    1:length(latent.by.year),
    function (idx) {
        reservoir.decay(
            latent.by.year[idx],
            6 * 30,
            365 * (10 - idx)
        )
    }
)

pdf("6mo_decay.pdf")
barplot(
    decay.6mo.hl / sum(decay.6mo.hl),
    names.arg=x.labels,
    main="Composition of reservoir by age (6mo half-life)",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir"
)
dev.off()

decay.1yr.hl <- sapply(
    1:length(latent.by.year),
    function (idx) {
        reservoir.decay(
            latent.by.year[idx],
            365,
            365 * (10 - idx)
        )
    }
)

pdf("1yr_decay.pdf")
barplot(
    decay.1yr.hl / sum(decay.1yr.hl),
    names.arg=x.labels,
    main="Composition of reservoir by age (1yr half-life)",
    xlab="Year prior to ART initiation",
    ylab="Proportion of reservoir"
)
dev.off()


