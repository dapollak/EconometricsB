library(lfe)
library(dplyr)
library(stargazer)
library(haven)

rm(list = ls())
options(scipen = 999)

#### Question 1 ####

N <- 1000000
mydata <- data.frame(matrix(, nrow = N, ncol = 0))

# This line generates an ability variable that takes 3 values:
mydata$ability <- floor(runif(N) * 3)

# This creates a female dummy which will only be used for part (l):
mydata$female <- as.numeric(rnorm(N) > 0)

# This line creates an indicator for living “near a college”.
# As written, the indicator will get the values 0 and 1 arbitrary.
# Note also that ability and near_col will not be correlated.
mydata$near_col <- as.numeric(rnorm(N) > 0)

# These two lines create a college attendance indicator.
# As constructed below, college attendance is somewhat correlated
# with both ability and being near a college.
# Make sure you exactly understand the construction
# of the “attend_college” indicator. Browse the data to verify.
mydata$attend_college <- 0
mydata$attend_college[mydata$ability == 2 |
        (mydata$ability == 1 & mydata$near_col == 1)] <- 1

# This last line generates white noise following a normal
# distribution with a standard deviation of 0.01 and a mean of zero;
mydata$white_noise <- rnorm(nrow(mydata), 0, 1)

# (b)
mydata$income <- 1 * mydata$ability + 2 * mydata$attend_college +
                            mydata$white_noise

# (c)
print(cor(mydata$ability, mydata$near_col))
print(cor(mydata$income, mydata$attend_college))

# (d)
reg1 <- lm(income ~ ability + attend_college, data = mydata)
print(summary(reg1))

# (e)
reg2 <- lm(income ~ attend_college, data = mydata)
print(summary(reg2))

# (f)
reg3_first_stage <- lm(attend_college ~ near_col, data = mydata)
print(summary(reg3_first_stage))
reg3 <- felm(income ~ 1 | 0 | (attend_college ~ near_col),
                        data = mydata)
print(summary(reg3))

# (g)
# (h)
    # compliers - mydata$ability == 1
    # always takers - mydata$ability == 2
    # never takers - mydata$ability == 0

# (i)
reg4 <- felm(income ~ 1 | 0 | (attend_college ~ near_col),
                            data = mydata[mydata$ability > 0, ])
print(summary(reg4))
print(stargazer(reg4, reg3, type = "text"))
# (j)
mydata$income_hetero <- 0
mydata <- mydata %>%
            group_by(ability) %>%
            mutate(income_hetero =
                    1 * ability +
                    (
                        ifelse(ability == 0, 1,
                        ifelse(ability == 1, 2, 10))
                    ) * attend_college + white_noise) %>%
                    ungroup()

# (k)
reg5 <- felm(income_hetero ~ 1 | 0 | (attend_college ~ near_col),
                            data = mydata[mydata$ability > 0, ])
stargazer(reg5, type = "text")

#### Question 2 ####
gelbach_data <- read_dta("data/gelbach.dta")
gelbach_data$bornearly <- ifelse(gelbach_data$quarter == 1, 1, 0)

# (c)
wald_estimator <- (mean(gelbach_data$hours[gelbach_data$bornearly == 1]) -
        mean(gelbach_data$hours[gelbach_data$bornearly == 0])) /
        (mean(gelbach_data$public[gelbach_data$bornearly == 1]) -
        mean(gelbach_data$public[gelbach_data$bornearly == 0]))
print(wald_estimator)

reg6_first_stage <- lm(
    public ~ bornearly,
    data = gelbach_data[gelbach_data$quarter < 2, ]
)
reg6 <- felm(hours ~ 1 | 0 | (public ~ bornearly),
                            data = gelbach_data[[gelbach_data$quarter < 2, ]])
stargazer(reg6, reg6_first_stage, type = "text")

# (e)
print(nrow(gelbach_data[gelbach_data$public == 1, ]))
print(nrow(gelbach_data[gelbach_data$bornearly == 1, ]))
print(nrow(gelbach_data[gelbach_data$public == 0, ]))
print(nrow(gelbach_data[gelbach_data$bornearly == 0, ]))

# (f)

# (i)
gelbach_data_youngest_is_5 <- gelbach_data[gelbach_data$youngest == 5, ]
gelbach_data_youngest_is_5$quarter <- ifelse(
    gelbach_data_youngest_is_5$quarter == 0,
    4,
    gelbach_data_youngest_is_5$quarter
)
balance_table_i <- gelbach_data_youngest_is_5 %>% select(
                    quarter, age, age2, grade, centcity, white,
                    num05, num612, num1317, othrlt18, othrge18
                    ) %>%
                    summarise(across(everything(), mean)) %>%
                    rbind(
                        gelbach_data_youngest_is_5 %>% select(
                            quarter, age, age2, grade, centcity, white,
                            num05, num612, num1317, othrlt18, othrge18
                        ) %>%
                        group_by(quarter) %>%
                        summarise(across(everything(), mean))
                    ) %>%
                    mutate(quarter = c(
                        "Full Sample", "74:II", "74:III", "74:IV", "75:I"
                    ))
balance_table_i <- t(balance_table_i)

# (j)
outcome_vars <- c("work79", "hours79", "weeksw79", "salary", "hours", "lfp")
results <- data.frame(matrix(nrow = 0, ncol = 3))

for (outcome in outcome_vars) {
    # outcome <- "work79"
    ols1_formula <- formula(sprintf("%s ~ public", outcome))
    ols1 <- lm(ols1_formula, data = gelbach_data_youngest_is_5)

    ols2_formula <- formula(sprintf("%s ~ public + age + age2 + grade + centcity + white + num05 + num612 + num1317 + othrlt18 + othrge18 | state", outcome))
    ols2 <- felm(ols2_formula, data = gelbach_data_youngest_is_5)

    iv_formula <- formula(
        sprintf("%s ~ age + age2 + grade + centcity + white + num05 + num612 + num1317 + othrlt18 + othrge18 | state | (public ~ quarter)",
        outcome)
    )
    iv <- felm(iv_formula, data = gelbach_data_youngest_is_5)

    results <- results %>% rbind(c(
        ols1$coefficients["public"],
        ols2$coefficients[1],
        iv$coefficients[11]
    ))
}

colnames(results) <- c("OLS1", "OLS2", "IV")
results
