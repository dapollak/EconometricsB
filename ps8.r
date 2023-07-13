library(lfe)
library(dplyr)
library(stargazer)
library(haven)
library(ggplot2)
library(binsreg)
library(rdrobust)

rm(list = ls())
options(scipen = 999)

#### Section 2 ####
data <- read_dta("data/adkw.dta")
data$VLBW <- ifelse(data$dbirwt < 1500, 1, 0)

## 2.2
reg1 <- lm(death1year ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500), data)

new_data <- data.frame(dbirwt = c(1450, 1550, 1499, 1501))
new_data$VLBW <- ifelse(new_data$dbirwt < 1500, 1, 0)
predict(reg1, new_data)

# 2.7
data$fitted_values <- as.numeric(predict(reg1, data))
g <- binsreg(y = data$death1year, x = data$dbirwt, nbins = 11)
g$bins_plot + geom_vline(xintercept = 1500, color = "red") + xlim(1350, 1650) + geom_line(data = data, aes(x = dbirwt, y = fitted_values))

#### Section 3 ####

# 3.1
g <- binsreg(y = data$death1year, x = data$dbirwt, nbins = 11)
g$bins_plot + 
    geom_vline(xintercept = 1500, color = "red") +
    xlab("grams") + ylab("mortality rate")

# 3.2
reg2 <- felm(death1year ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500) + gestat | 0 | 0 | dbirwt, data)
print(stargazer(reg2, type = "text"))

# 3.7 (i) and (ii)
ggplot(data, aes(x=dbirwt)) + geom_histogram()
data <- data %>%
        group_by(dbirwt) %>%
        mutate(count = n()) %>%
        ungroup()

data_test <- data %>%
            select(c("count", "dbirwt", "VLBW")) %>%
            group_by(dbirwt) %>%
            summarise(
                count = min(count),
                dbirwt = min(dbirwt),
                VLBW = min(VLBW)
            )

reg3 <- lm(count ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500), data_test)
stargazer(reg3, type = "text")

# 3.8
g <- binsreg(y = data$gestat, x = data$dbirwt, nbins = 11)
g$bins_plot + geom_vline(xintercept = 1500, color = "red") +
    xlab("grams") + ylab("gestational age")

# 3.9
reg4 <- lm(death1year ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500) + VLBW:I((dbirwt - 1500)^2) + I(1-VLBW):I((dbirwt - 1500)^2), data)
reg5 <- lm(death1year ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500) + VLBW:I((dbirwt - 1500)^2) + I(1-VLBW):I((dbirwt - 1500)^2) + VLBW:I((dbirwt - 1500)^3) + I(1-VLBW):I((dbirwt - 1500)^3), data)
reg6 <- lm(death1year ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500) + VLBW:I((dbirwt - 1500)^2) + I(1-VLBW):I((dbirwt - 1500)^2) + VLBW:I((dbirwt - 1500)^3) + I(1-VLBW):I((dbirwt - 1500)^3) + VLBW:I((dbirwt - 1500)^4) + I(1-VLBW):I((dbirwt - 1500)^4), data)
stargazer(reg4, reg5, reg6, type = "text")

# plot quad
data$fitted_values_quad <- as.numeric(predict(reg4, data))
g <- binsreg(y = data$death1year, x = data$dbirwt, nbins = 11)
g$bins_plot + geom_vline(xintercept = 1500, color = "red") + xlim(1350, 1650) + geom_line(data = data, aes(x = dbirwt, y = fitted_values_quad))

# 3.10
reg7 <- rdrobust(y = data$death1year, x = data$dbirwt, c = 1500, p=2, q=3, vce = "hc0", cluster=data$dbirwt, bwselect = "mserd", kernel = "Triangular")
reg8 <- rdrobust(y = data$death1year, x = data$dbirwt, c = 1500, p=3, q=4, vce = "hc0", cluster=data$dbirwt, bwselect = "mserd", kernel = "Triangular")
reg9 <- rdrobust(y = data$death1year, x = data$dbirwt, c = 1500, p=4, q=5, vce = "hc0", cluster=data$dbirwt, bwselect = "mserd", kernel = "Triangular")

# 3.11
estimates_bw_df <- data.frame(bandwidth = seq(85, 5, -5), estimates = 0, ci_lower = 0, ci_upper = 0)
for (bw in estimates_bw_df$bandwidth) {
    reg_bw <- rdrobust(y = data$death1year, x = data$dbirwt, c = 1500, p=1, q=2, h = bw, vce = "hc0", cluster=data$dbirwt)
    estimates_bw_df[estimates_bw_df$bandwidth == bw,]$estimates <- reg_bw$Estimate[1]
    estimates_bw_df[estimates_bw_df$bandwidth == bw,]$ci_lower <- reg_bw$ci[1]
    estimates_bw_df[estimates_bw_df$bandwidth == bw,]$ci_upper <- reg_bw$ci[4]
}

for (bw in estimates_bw_df$bandwidth) {
    bw_data <- subset(data, (dbirwt > (1500 - bw/2)) & (dbirwt <= (1500 + bw/2)))
    reg_bw <- felm(death1year ~ VLBW + VLBW:I(dbirwt - 1500) + I(1-VLBW):I(dbirwt - 1500) + gestat | 0 | 0 | dbirwt, bw_data)
    estimates_bw_df[estimates_bw_df$bandwidth == bw,]$estimates <- reg_bw$coefficients[2]
    estimates_bw_df[estimates_bw_df$bandwidth == bw,]$ci_lower <- confint(reg_bw, "VLBW")[1]
    estimates_bw_df[estimates_bw_df$bandwidth == bw,]$ci_upper <- confint(reg_bw, "VLBW")[2]
}

ggplot(estimates_bw_df, aes(x = bandwidth, y = estimates)) + geom_point() +
geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), colour="blue", width=.1)
