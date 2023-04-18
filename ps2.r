library(dplyr)
library(plm)
library(haven)
options(scipen=999)

simulate_data <- function(
    individuals, periods, alpha, beta, rho_u, rho_tau, tau_vec, gamma = 0
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
        mutate(y = alpha + gamma * year + beta * x + c + u)

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
    critical_value,
    gamma = 0
) {
    tau <-  runif(individuals, 2, periods + 1)
    results <- data.frame(matrix(NA, nrow = repeats, ncol = 7))
    colnames(results) <- c(
        "fe_no_corr", "fe_corr", "fd_no_corr", "fd_corr", "polled_beta", "fe_beta", "fd_beta"
    )

    for (r in 1:repeats) {
        dat <- simulate_data(
            individuals, periods, alpha, beta, rho_u, rho_tau, tau, gamma
        )

        ## Pooled OLS
        results$polled_beta[r] <- lm(y ~ x, data = dat)$coefficients["x"]

        pdat <- pdata.frame(dat, index = c("id", "year"))

        ## First Difference
        fd_reg1 <- plm(y ~ x, data = pdat, model = "fd")
        results$fd_beta[r] <- fd_reg1$coefficients["x"]

        # No Clustered SE
        results$fd_no_corr[r] <-
            abs(summary(fd_reg1)$coefficients[6]) >= critical_value

        # Clustered SE
        rob_vcov2 <- vcovHC(fd_reg1, type = "sss", cluster = "group")
        results$fd_corr[r] <-
            abs(summary(fd_reg1, rob_vcov2)$coefficients[6]) >= critical_value

        ## Fixed Effects
        fe_reg1 <- plm(
            y ~ x, data = pdat, model = "within", effect = "individual"
        )
        results$fe_beta[r] <- fe_reg1$coefficients["x"]

        # No Clustered SE
        results$fe_no_corr[r] <-
            abs(summary(fe_reg1)$coefficients[3]) >= critical_value

        # Clustered SE
        rob_vcov1 <- vcovHC(fe_reg1, type = "sss", cluster = "group")
        results$fe_corr[r] <-
            abs(summary(fe_reg1, rob_vcov1)$coefficients[3]) >= critical_value
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

# (2 + 3)
r1 <- monte_carlo(250, 20, 5, 0, 0, 0.02, 200, 1.96)
print(c(
    (rmse_bias(0, r1$polled_beta))^2,
    rmse_variation(0, r1$polled_beta),
    rmse_calc(0, r1$polled_beta)
))
r2 <- monte_carlo(250, 20, 5, 0, 0.9, 0.02, 200, 1.96)


# (4)
r3 <- monte_carlo(250, 20, 5, 0.5, 0, 0.02, 200, 1.96)
print(c(
    (rmse_bias(0.5, r3$fe_beta))^2,
    rmse_variation(0.5, r3$fe_beta),
    rmse_calc(0.5, r3$fe_beta)
))
print(c(
    (rmse_bias(0.5, r3$fd_beta))^2,
    rmse_variation(0.5, r3$fd_beta),
    rmse_calc(0.5, r3$fd_beta)
))
r4 <- monte_carlo(250, 20, 5, 0.5, 0.9, 0.02, 200, 1.96)
print(c(
    (rmse_bias(0.5, r4$fe_beta))^2,
    rmse_variation(0.5, r4$fe_beta),
    rmse_calc(0.5, r4$fe_beta)
))
print(c(
    (rmse_bias(0.5, r4$fd_beta))^2,
    rmse_variation(0.5, r4$fd_beta),
    rmse_calc(0.5, r4$fd_beta)
))

# Part B, (3)
r5 <- monte_carlo(250, 20, 5, 0.5, 0, 0.02, 200, 1.96, 0.05)
r6 <- monte_carlo(250, 20, 5, 0.5, 0.9, 0.02, 200, 1.96, 0.05)


####### Question 2

immg_data <- read_dta("data/immigration8994.dta")

reg_formula <- lnhwage ~ imm3 + sch + exp +
                expsq + married + nonjew + for_born + timeisr +
                year90  + year91 + year92 + year93 + year94

immg_pdata <- pdata.frame(immg_data, index = c("occ2", "year"))
immg_pooled_reg <- lm(reg_formula, data = immg_data) # beta = -0.2083 (0.182)

immg_fe_reg <- plm(reg_formula, data = immg_pdata,
                    model = "within", effect = "individual") # beta = -0.29550569 (0.10204864) cl - 0.10665030

immg_fd_reg <- plm(reg_formula, data = immg_pdata, model = "fd") # beta = 0.046 (0.18408314) cl - 0.28318235

immg_fe_vcov <- vcovHC(immg_fe_reg, type = "sss", cluster = "group")
immg_fd_vcov <- vcovHC(immg_fd_reg, type = "sss", cluster = "group")

## Serial correlation test
pwfdtest(immg_fd_reg, h0 = "fd")