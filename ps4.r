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

# d <- simulate_data(50, 8, 1, 0.25, 0.25, 0.25)
r1 <- monte_carlo(50, 8, 1, 0.25, 0.25, 0.25, 200)
rmse1 <- calc_rmse(0.25, r1$with_controls_gamma)
rmse2 <- calc_rmse(0.25, r1$no_controls_gamma)

# 2(a)
r2 <- monte_carlo(50, 8, 0, 0.25, 0.25, 0.25, 200,
                    throw_lines = floor(0.25 * 50))
rmse3 <- calc_rmse(0.25, r2$no_controls_gamma)

# 2(b)
# d <- simulate_data(50, 8, 0, 0.25, 0.25, 0.25, throw_per_of_control = 0.25)
r3 <- monte_carlo(50, 8, 0, 0.25, 0.25, 0.25, 200, throw_per_of_control = 0.25)
rmse4 <- calc_rmse(0.25, r3$no_controls_gamma)

# 3(a)
r4 <- monte_carlo(50, 8, 0, 0.25, 0.25, 0.25, 200,
                throw_per_of_control = 0.25, throw_per_of_treat = 0.5)
rmse5 <- calc_rmse(0.25, r4$no_controls_gamma)


function_a <- function(a, b, c) {
    print("hello")
    data_in_function <- c(1, 2, 4)
    print(data_in_function)
    print(b)
    print("bye")
}

function_a(1, 2, 30)
function_a(10, 25, 34)
