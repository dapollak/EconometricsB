library(dplyr)
library(plm)
library(haven)
options(scipen = 999)

simulate_data <- function(
    individuals,
    alpha,
    beta,
    gamma,
    sig_u,
    sig_x,
    throw_lines = 0,
    throw_per_of_control = 0,
    throw_per_of_treat = 0
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

    if (throw_lines > 0) {
        dat <- dat[-sample(1:individuals, throw_lines), ]
    }

    if (throw_per_of_control > 0) {
        per <- quantile(dat$y, throw_per_of_control)

        if (throw_per_of_treat) {
            rows_to_drop <- dat[dat$y < per & dat$t == 1, ]
            rows_to_drop <- rows_to_drop[sample(
                nrow(rows_to_drop),
                floor(throw_per_of_treat * nrow(rows_to_drop))
            ), ]
            if (nrow(rows_to_drop)) {
                row_indices_to_drop <- as.numeric(rownames(rows_to_drop))
                dat <- dat[-row_indices_to_drop, ]
            }
        } else {
            rows_to_drop <- dat[dat$t == 0 & dat$y < per, ]
            row_indices_to_drop <- as.numeric(rownames(rows_to_drop))
            dat <- dat[-row_indices_to_drop, ]
        }
    }

    return(dat)
}

monte_carlo <- function(
    individuals,
    alpha,
    beta,
    gamma,
    sig_u,
    sig_x,
    repeats,
    throw_lines = 0,
    throw_per_of_control = 0,
    throw_per_of_treat = 0
) {
    results <- data.frame(matrix(NA, nrow = repeats, ncol = 4))
    colnames(results) <- c(
        "no_controls_gamma", "with_controls_gamma",
        "no_controls_reject", "with_controls_reject"
    )

    for (r in 1:repeats) {
        dat <- simulate_data(
            individuals, alpha, beta, gamma, sig_u,
            sig_x, throw_lines, throw_per_of_control, throw_per_of_treat
        )

        r1 <- lm(y ~ x + t, data = dat)
        r2 <- lm(y ~ t, data = dat)

        results$with_controls_gamma[r] <- r1$coefficients["t"]
        results$no_controls_gamma[r] <- r2$coefficients["t"]

        # t-values
        results$no_controls_reject <- summary(r1)$coefficients[, 3]["t"] > 1.96
        results$with_controls_reject <-
                summary(r2)$coefficients[, 3]["t"] > 1.96
    }

    return(results)
}

calc_rmse <- function(parameter, estimator_vector) {
    results <- data.frame(matrix(NA, nrow = 1, ncol = 3))
    colnames(results) <- c(
        "rmse", "bias", "variation"
    )

    results$rmse <- (1 / length(estimator_vector)) *
            sum((estimator_vector - parameter)^2)
    results$bias <- (mean(estimator_vector) - parameter)
    results$variation <- (1 / length(estimator_vector)) *
            (sum((mean(estimator_vector) - estimator_vector)^2))

    return(results)
}

build_residual_diff <- function(N, T, fd_residuals) {
    result <- data.frame(matrix(NA, nrow = N * (T - 2), ncol = 2))
    colnames(result) <- c("du", "dup")
    for (i in 1:N) {
        for (j in 1:(T - 2)) {
            result[(1 + (i - 1) * (T - 2)) + (j - 1), ] <-
                    c(
                        fd_residuals[(1 + (i - 1) * (T - 1)) + j],
                        fd_residuals[(1 + (i - 1) * (T - 1)) + (j - 1)]
                    )
        }
    }

    return(result)
}

build_residual_diff2 <- function(N, T, fd_residuals) {
    result <- data.frame(matrix(NA, nrow = N * (T - 1), ncol = 2))
    colnames(result) <- c("u", "up")
    for (i in 1:N) {
        for (j in 1:(T - 1)) {
            result[(1 + (i - 1) * (T - 2)) + (j - 1), ] <-
                    c(
                        fd_residuals[(1 + (i - 1) * T) + j],
                        fd_residuals[(1 + (i - 1) * T) + (j - 1)]
                    )
        }
    }

    return(result)
}

#### Question 2 ####
pre_data <- read_dta("data/preexperiment.dta")
pre_pdat <- pdata.frame(pre_data, index = c("id", "t"))
fd_reg1 <- plm(y ~ x, data = pre_pdat, model = "fd")
beta1 <- fd_reg1$coefficients["x"]

rr <- build_residual_diff(1431, 4, fd_reg1$residuals)
res_reg <- lm(du ~ dup - 1, data = rr)
rho <- res_reg$coefficients["dup"]

## simple way no FE

reg2 <- lm(y ~ x, data = pre_data)
beta0 <- reg2$coefficients[1]
beta1 <- reg2$coefficients["x"]
rr <- build_residual_diff2(1431, 4, reg2$residuals)
res_reg <- lm(u ~ up - 1, data = rr)
rho <- res_reg$coefficients["up"]
var_epsilon <- var(res_reg$residuals)
var_x <- var(pre_data["x"])