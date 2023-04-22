library(dplyr)
library(plm)
library(haven)
options(scipen = 999)

simulate_data <- function(
    individuals, alpha, beta, gamma, sig_u, sig_x
) {
    # data frame
    dat <- expand.grid(1:individuals)
    colnames(dat) <- c("id")

    # Generate c_i, x_it and u_it
    dat <- dat %>%
        arrange(id) %>%
        mutate(
            x = rnorm(individuals, 0, sig_x),
            t = sample(c(0, 1), individuals, replace = TRUE),
            u = rnorm(individuals, 0, sig_u)
        )

    # Generate y according to the DGP
    dat <- dat %>%
        mutate(
            y = alpha + (beta * x) + (gamma * t) + u,
        )

    return(dat)
}

monte_carlo <- function(
    individuals,
    alpha,
    beta,
    gamma,
    sig_u,
    sig_x,
    repeats
) {
    results <- data.frame(matrix(NA, nrow = repeats, ncol = 2))
    colnames(results) <- c(
        "no_controls_beta", "with_controls_beta"
    )

    for (r in 1:repeats) {
        dat <- simulate_data(
            individuals, alpha, beta, gamma, sig_u, sig_x
        )

        results$with_controls_beta[r] <- 
                    lm(y ~ x + t, data = dat)$coefficients["t"]
        results$no_controls_beta[r] <- lm(y ~ t, data = dat)$coefficients["t"]
    }

    return(results)
}

rmse_calc <- function(beta, beta_vector) {
    R <- length(beta_vector)
    result <- (1 / R) * sum((beta_vector - beta)^2)
    return(result)
}

rmse_bias <- function(beta, beta_vector) {
    return(mean(beta_vector) - beta)
}

rmse_variation <- function(beta, beta_vector) {
    return((1 / length(beta_vector)) *
            (sum((mean(beta_vector) - beta_vector)^2)))
}

# d <- simulate_data(50, 8, 1, 0.25, 0.25, 0.25)
r <- monte_carlo(50, 8, 1, 0.25, 0.25, 0.25, 200)