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

#### A2
monte_carlo_pooled <- function(
    individuals,
    periods,
    alpha,
    beta,
    rho_u,
    rho_tau,
    repeats
) {
    tau <-  runif(individuals, 2, periods + 1)
    betas <- rep(0, repeats)
    for (r in 1:repeats) {
        dat <- simulate_data(
            individuals, periods, alpha, beta, rho_u, rho_tau, tau
        )

        betas[r] <- lm(y ~ x, data = dat)$coefficients["x"]
    }

    return(betas)
}
# betas <- monte_carlo_pooled(250, 20, 5, 0, 0, 0.02, 200)

#### A3
monte_carlo_fe <- function(
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
    rejection_table <- data.frame(matrix(NA, nrow = repeats, ncol = 4))
    colnames(rejection_table) <- c(
        "fe_no_corr", "fe_corr", "fd_no_corr", "fd_corr"
    )

    for (r in 1:repeats) {
        dat <- simulate_data(
            individuals, periods, alpha, beta, rho_u, rho_tau, tau
        )

        pdat <- pdata.frame(dat, index = c("id", "year"))

        ## Fixed Effects
        fe_reg1 <- plm(
            y ~ x, dat = pdat, model = "within", effect = "individual"
        )

        # No Clustered SE
        rejection_table$fe_no_corr[r] <-
            abs(summary(fe_reg1)$coefficients[3]) >= critical_value

        # Clustered SE
        rob_vcov1 <- vcovHC(fe_reg1, type = "sss", cluster = "group")
        rejection_table$fe_corr[r] <-
            abs(summary(fe_reg1, rob_vcov1)$coefficients[3]) >= critical_value

        ## First Difference
        fd_reg1 <- plm(y ~ x, dat = pdat, model = "fd")

        # No Clustered SE
        rejection_table$fd_no_corr[r] <-
            abs(summary(fd_reg1)$coefficients[3]) >= critical_value

        # Clustered SE
        rob_vcov2 <- vcovHC(fd_reg1, type = "sss", cluster = "group")
        rejection_table$fd_corr[r] <-
            abs(summary(fd_reg1, rob_vcov2)$coefficients[3]) >= critical_value
    }

    return(rejection_table)
}

rejections <- monte_carlo_fe(250, 20, 5, 0, 0, 0.02, 200, 1.96)