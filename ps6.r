library(dplyr)
library(plm)
library(lfe)
library(haven)
library(ggplot2)
library(did)
library(rstatix)
library(ggpubr)
library(tidyr)
library(MatchIt)
options(scipen = 999)

#### Question 4


### Section 1
nsw_dw <- read_dta("data/nsw_dw.dta")

balance_table <- nsw_dw %>%
                select(
                    treat, age, education, black, hispanic,
                    married, nodegree, re74, re75, u74
                ) %>%
                pivot_longer(
                    -treat,
                    names_to = "variable",
                    values_to = "value"
                ) %>%
                group_by(variable) %>%
                t_test(value ~ treat, detailed = T)

balance_reg <- lm(education + age ~ treat, data = nsw_dw)

print(mean(nsw_dw[nsw_dw$treat == 1, ]$re78) -
      mean(nsw_dw[nsw_dw$treat == 0, ]$re78))

sec1_reg <- lm(re78 ~ treat + age + education + black + hispanic +
                married + nodegree + re74 + re75 + u74, data = nsw_dw)
summary(sec1_reg)
### Section 2

nsw_psid_dw <- read_dta("data/nsw_psid_dw.dta")
balance_table2 <- nsw_psid_dw %>%
                select(
                    treat, age, education, black, hispanic,
                    married, nodegree, re74, re75, u74
                ) %>%
                pivot_longer(
                    -treat,
                    names_to = "variable",
                    values_to = "value"
                ) %>%
                group_by(variable) %>%
                t_test(value ~ treat, detailed = T)

# As calculated in Dehejia, Rajeev and S. Wahba (2002)
print(mean(nsw_psid_dw[nsw_psid_dw$treat == 1, ]$re78) -
      mean(nsw_psid_dw[nsw_psid_dw$treat == 0, ]$re78))

sec2_reg <- lm(re78 ~ treat + age + education + black + hispanic +
                married + nodegree + re74 + re75 + u74, data = nsw_psid_dw)
summary(sec2_reg)

### Section 3

# calc propensity score
nsw_psid_dw <- nsw_psid_dw %>%
            mutate(
                age_sq = age^2,
                education_sq = education^2,
                re74_sq = re74^2,
                re75_sq = re75^2
            )

propensity_score_lm <- glm(treat ~ age + age_sq + education + education_sq +
                        black + hispanic + married + nodegree +
                        re74 + re74_sq + re75 + re75_sq + u74 +
                        black:u74, data = nsw_psid_dw,
                        family = binomial("logit"))
summary(propensity_score_lm)
nsw_psid_dw$pscore <- predict(propensity_score_lm, type = "response")

# 8
ggplot(nsw_psid_dw, aes(x = pscore, fill = as.factor(treat))) +
    geom_histogram()

# 9
min_max_tab <- nsw_psid_dw %>%
                group_by(treat) %>%
                summarize(minima = min(pscore), maxima = max(pscore))
nsw_psid_dw <- nsw_psid_dw %>%
        filter(pscore < min(min_max_tab$maxima) &
                    pscore > max(min_max_tab$minima))

# 10

balance_table3 <- nsw_psid_dw %>%
                select(
                    treat, age, education, black, hispanic,
                    married, nodegree, re74, re75, u74
                ) %>%
                pivot_longer(
                    -treat,
                    names_to = "variable",
                    values_to = "value"
                ) %>%
                group_by(variable) %>%
                t_test(value ~ treat, detailed = T)

balance_reg3 <- lm(education + age ~ treat, data = nsw_psid_dw)

print(mean(nsw_psid_dw[nsw_psid_dw$treat == 1, ]$re78) -
      mean(nsw_psid_dw[nsw_psid_dw$treat == 0, ]$re78))

sec3_reg <- lm(re78 ~ treat + age + education + black + hispanic +
                married + nodegree + re74 + re75 + u74, data = nsw_dw)
summary(sec1_reg)

# 11 - NN Matching
for (i in c(1, 5)) {
    nn_match <- matchit(treat ~ age + age_sq + education + education_sq +
                            black + hispanic + married + nodegree +
                            re74 + re74_sq + re75 + re75_sq + u74 +
                            black:u74,
                            data = nsw_psid_dw,
                            method = "nearest",
                            distance = "mahalanobis",
                            estimand = "ATT",
                            replace = T,
                            caliper = NULL,
                            exact = NULL,
                            ratio = i)

    nn_match_data <- match.data(nn_match)
    nn_match_reg <- lm(
        re78 ~ treat,
        data = nn_match_data,
        weights = nn_match_data$weights
    )

    print(summary(nn_match_reg)$coefficients["treat", 1])
}

# 13
prop_match <- matchit(treat ~ age + age_sq + education + education_sq +
                            black + hispanic + married + nodegree +
                            re74 + re74_sq + re75 + re75_sq + u74 +
                            black:u74,
                            data = nsw_psid_dw,
                            method = "nearest",
                            distance = "glm",
                            link = "logit",
                            estimand = "ATT",
                            replace = T,
                            caliper = NULL,
                            exact = NULL,
                            ratio = 5)
prop_match_data <- match.data(prop_match)
prop_match_reg <- lm(
    re78 ~ treat,
    data = prop_match_data,
    weights = prop_match_data$weights
)

print(summary(prop_match_reg)$coefficients["treat", 1])
