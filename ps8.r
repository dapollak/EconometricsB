library(lfe)
library(dplyr)
library(stargazer)
library(haven)
library(ggplot2)
library(binsreg)

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
g$bins_plot + geom_vline(xintercept = 1500, color = "red") + xlab("X") + ylab("Y") + xlim(1350, 1650)
ggplot(data = data, aes(x = dbirwt, y = fitted_values)) + geom_line() + xlab("grams") + ylab("mortality rate") + xlim(0, 3000)

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
