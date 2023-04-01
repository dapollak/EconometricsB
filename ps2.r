library(dplyr)
library(plm)

simulate_data <- function(
    individuals, periods, alpha, beta, rho_u, rho_tau, tau_vec
) {
    # Panel data frame
    dat <- expand.grid(1:individuals, 1:periods)
    colnames(dat) <- c("id", "year")

    # Function to simulate u_it
    simulate_error <- function(n_periods, rho = 0, sd = 0.1) {
        u <- rep(0, n_periods + 1)
        for (t in 1:n_periods) {
            u[t + 1] <- rho * u[t] + rnorm(1, 0, sd)
        }
        return(u[-1])
    }

    # Generate c_i, x_it and u_it
    dat <- dat %>%
        arrange(id, year) %>%
        group_by(id) %>%
        mutate(
            tau = tau_vec[cur_group_id()],
            c = rho_tau * tau + rnorm(1, 0, 0.1),
            x = ifelse(year >= floor(tau), 1, 0),
            u = simulate_error(n_periods = n(), rho = rho_u)
        ) %>%
        ungroup()
    # Generate y according to the DGP
    dat <- dat %>%
        mutate(y = alpha + beta * x + c + u)

    return(dat)
}

monte_carlo <- function(
    individuals,
    periods,
    alpha,
    beta,
    rho_u,
    rho_tau,
    repeats,
    critical_value
) {
    tau <-  runif(individuals, 2, periods + 1)
    results <- data.frame(matrix(NA, nrow = repeats, ncol = 7))
    colnames(results) <- c(
        "fe_no_corr", "fe_corr", "fd_no_corr", "fd_corr", "polled_beta", "fe_beta", "fd_beta"
    )

    for (r in 1:repeats) {
        dat <- simulate_data(
            individuals, periods, alpha, beta, rho_u, rho_tau, tau
        )

        ## Pooled OLS
        results$polled_beta[r] <- lm(y ~ x, data = dat)$coefficients["x"]

        pdat <- pdata.frame(dat, index = c("id", "year"))

        ## Fixed Effects
        fe_reg1 <- plm(
            y ~ x, data = dat, model = "within", effect = "individual"
        )
        results$fe_beta[r] <- fe_reg1$coefficients["x"]

        # No Clustered SE
        results$fe_no_corr[r] <-
            abs(summary(fe_reg1)$coefficients[3]) >= critical_value

        # Clustered SE
        rob_vcov1 <- vcovHC(fe_reg1, type = "sss", cluster = "group")
        results$fe_corr[r] <-
            abs(summary(fe_reg1, rob_vcov1)$coefficients[3]) >= critical_value

        ## First Difference
        fd_reg1 <- plm(y ~ x, data = pdat, model = "fd")
        results$fd_beta[r] <- fd_reg1$coefficients["x"]

        # No Clustered SE
        results$fd_no_corr[r] <-
            abs(summary(fd_reg1)$coefficients[3]) >= critical_value

        # Clustered SE
        rob_vcov2 <- vcovHC(fd_reg1, type = "sss", cluster = "group")
        results$fd_corr[r] <-
            abs(summary(fd_reg1, rob_vcov2)$coefficients[3]) >= critical_value
    }

    return(results)
}

r <- monte_carlo(250, 20, 5, 0, 0, 0.02, 200, 1.96)