library(lfe)
library(dplyr)
library(stargazer)

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
