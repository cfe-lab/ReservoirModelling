# In this file, we break out and clean up some of the code in ode.r and turn it
# into a library.

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


reservoir.decay <- function(size, half.life, decay.time) {
    # Parameters:
    # size is the number of particles in question.
    # half.life is the half-life in days.
    # decay.time is the amount of time the particles have aged.
    return(size * 2^(-(decay.time / half.life)))
}


bin.helper <- function(
    ode.out.df,
    bin.size=365,
    grid.size=0.001
) {
    # Helper to assist with binning the latent reservoir into time periods.
    bin.indices <- seq(nrow(ode.out.df), 1, by=-(bin.size * (1 / grid.size)))
    if (bin.indices[length(bin.indices)] != 1) {
        bin.indices <- c(bin.indices, 1)
    }
    bin.indices <- rev(bin.indices)

    latent.by.bin <- (
        ode.out.df$L.P[bin.indices[2:length(bin.indices)]] - 
        ode.out.df$L.P[bin.indices[1:length(bin.indices)-1]]
    )
    return(latent.by.bin)
}


undecayed.reservoir.distribution <- function(
    days.pre.therapy,
    bin.size=365, 
    grid.size=0.001, 
    maxsteps=25000
) {
    # Get an undecayed "distribution" of the reservoir composition by age.
    # days.pre.therapy is the amount of time the infection was unchecked.
    # bin.size is the size of each "time bin" in the distribution (in "grid points").
    # grid.size is the step size used in the numerical solution of the ODEs defining 
    #     the dynamics of the untreated infection.  (To make things simpler, this should
    #     be a number such that 1 / grid.size is an integer.)
    times <- seq(0, days.pre.therapy, by=grid.size)

    out <- ode(
        y=state,
        times=times,
        func=derivatives,
        parms=parameters,
        maxsteps=maxsteps
    )

    out.df <- as.data.frame(out)
    latent.by.bin <- bin.helper(out.df, bin.size=bin.size, grid.size=grid.size)
    
    return(
        list(
            solution=out.df,
            bin.freqs=latent.by.bin,
            bin.dist.no.decay=latent.by.bin / sum(latent.by.bin)
        )
    )
}


undecayed.reservoir.distribution.given.vl <- function(
    vl.times,  # given as days relative to (estimated) infection date
    vl.measurements,
    infection.date,  # date object
    art.initiation,  # date object
    bin.size=365, 
    grid.size=0.001, 
    maxsteps=25000
) {
    virus <- approxfun(vl.times, vl.measurements, rule=2)

    state.vl.known <- state[1:6]

    derivatives.vl.known <- function(t, state, parameters) {
        with(
            as.list(c(state, parameters)), 
            {
                V <- virus(t)

                dS <- alpha.S - delta.S * S - beta * S * V
                dA.P <- (1 - lambda) * tau * beta * S * V - delta.I * A.P - kappa * A.P * E
                dA.U <- (1 - lambda) * (1 - tau) * beta * S * V - delta.I * A.U - kappa * A.U * E
                dL.P <- lambda * tau * beta * S * V
                dL.U <- lambda * (1 - tau) * beta * S * V
                dE <- alpha.E + omega * (A.P + A.U) * E / (E + E.50) - delta.E * E

                return(
                    list(
                        c(dS, dA.P, dA.U, dL.P, dL.U, dE),
                        V=V
                    )
                )
            }
        )
    }

    max.time <- as.numeric(art.initiation - infection.date, units="days")
    times <- seq(0, max.time, by=grid.size)

    run.time <- system.time(
        vl.known <- ode(
            y=state.vl.known,
            times=times,
            func=derivatives.vl.known,
            parms=parameters,
            maxsteps=25000
        )
    )
    cat("Time elapsed:\n")
    print(run.time)

    out.df <- as.data.frame(vl.known)
    latent.by.bin <- bin.helper(out.df, bin.size=bin.size, grid.size=grid.size)
    
    return(
        list(
            solution=out.df,
            bin.freqs=latent.by.bin,
            bin.dist.no.decay=latent.by.bin / sum(latent.by.bin)
        )
    )
}


decay.distribution <- function(
    latent.by.bin,
    half.life,
    bin.size
) {
    # Decay the given undecayed distribution given a reservoir half-life.
    # half.life is the half life (in days) of reservoir particles.
    # If there are, e.g., 10 bins, they look like
    # [10-9 years prior, 9-8 years prior, ..., 1-0 years prior]
    # and we decay each one by 
    # [9 years, 8 years, ,.., 0 years]
    # If the earliest one is not a full year we still decay it from the
    # end of the time period.
    decayed <- sapply(
        1:length(latent.by.bin),
        function (idx) {
            reservoir.decay(
                latent.by.bin[idx],
                half.life,
                bin.size * (length(latent.by.bin) - idx)
            )
        }
    )
    
    return(
        list(
            bin.freqs=decayed,
            bin.dist=decayed / sum(decayed)
        )
    )
}


get.bin.helper <- function(seq.age, bin.size) {
    floor(seq.age / bin.size)
}


log.likelihood.no.factorial <- function(seq.ages, bin.dist, bin.size) {
    likelihoods <- sapply(
        seq.ages,
        function (x) bin.dist[max(1, length(bin.dist) - get.bin.helper(x, bin.size))]
    )
    return(sum(log(likelihoods)))
}


fisher.information <- function(decay, bin.freqs) {
    # Compute the Fisher information for the given reservoir distribution and decay rate.
    # Note that the unit for decay is the bin size of bin.freqs, *not* necessarily days.
    n <- length(bin.freqs)

    h <- function(decay) {
        # This is the normalizing factor for the density.
        summands <- sapply(
            seq(1, n),
            function (j) {
                bin.freqs[j] * 2^(-(n - j) / decay)
            }
        )
        return(sum(summands))
    }

    h.prime <- function(theta) { 
        # This is the derivative of the above.
        summands <- sapply(
            seq(1, n),
            function (j) {
                bin.freqs[j] * log(2) * (n - j) * 2^(-(n-j) / theta) / theta^2
            }
        )
        return(sum(summands))
    }

    f <- function(x, theta) {
        numerator <- bin.freqs[x] * 2^(-(n-x) / theta)
        return(numerator / h(theta))
    }

    ll.prime <- function(x, theta) {
        first.term <- log(2) * (n - x) / theta^2
        second.term <- h.prime(theta) / h(theta)
        return(first.term - second.term)
    }

    # The Fisher information is \E_{\theta}[ (\frac{\partial}{\partial\theta} log f(X|\theta))^2 ],
    # where X is distributed as f(\cdot|\theta).
    summands <- sapply(
        seq(1, n),
        function (j) f(j, decay) * ll.prime(j, decay)^2
    )
    return(sum(summands))
}

