library(dplyr)
library(plm)

simulate_data <- function(
    individuals, periods, alpha, beta, rho_u, rho_tau, tau_vec = NULL
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
            tau = ifelse(is.null(tau_vec), runif(1, 2, periods + 1), tau_vec),
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
dat <- simulate_data(250, 20, 5, 0, 0, 0.02)

#### A3
pdat <- pdata.frame(dat, index = c("id", "year"))

## Fixed Effects
fe_reg1 <- plm(y ~ x, dat = pdat, model = "within", effect = "individual")

# No Clustered SE
summary(fe_reg1)

# Clustered SE
rob_vcov1 <- vcovHC(fe_reg1, type = "sss", cluster = "group")
summary(fe_reg1, rob_vcov1)

## First Difference
fd_reg1 <- plm(y ~ x, dat = pdat, model = "fd")
# No Clustered SE
summary(fd_reg1)

# Clustered SE
rob_vcov2 <- vcovHC(fd_reg1, type = "sss", cluster = "group")
summary(fd_reg1, rob_vcov2)