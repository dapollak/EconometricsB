library(dplyr)
library(plm)
library(lfe)
library(haven)
library(ggplot2)
library(did)
options(scipen = 999)

#### Question 1

indo_dat <- read_dta("data/Indo_Schooling.dta")
## (A)
summary(indo_dat)

# size of cohorts
print(
    indo_dat %>%
    group_by(birth_year) %>%
    summarise(num_of_people = n())
)

# size of regions
print(
    indo_dat %>%
    group_by(birth_region) %>%
    summarise(
        num_of_people = n(),
        num_of_schools = indo_dat[indo_dat$birth_region == 
                            birth_region, ]$num_schools[1],
        intense = indo_dat[indo_dat$birth_region ==
                            birth_region, ]$program_intensity[1]
    )
)
## (B)

reg1 <- lm(log_wage ~ education, data = indo_dat)

## (C)
intense_levels <- indo_dat[c("program_intensity", "num_schools", "education")]

# (C-a)
intense_levels <- intense_levels %>%
    group_by(program_intensity) %>%
    summarise(
        num_schools_avg = mean(num_schools),
        education_avg = mean(education),
        num_of_people = n()
    )
print(intense_levels)

# (C-b)
print(
    mean(indo_dat[indo_dat$program_intensity == 1, ]$education) -
    mean(indo_dat[indo_dat$program_intensity == 0, ]$education)
)

summary(lm(education ~ program_intensity, data = indo_dat))

# (C-c)
indo_dat_treated <- indo_dat %>%
    mutate(
        treated = ifelse(birth_year < 62, 0, 1)
    )

print(
    mean(indo_dat_treated[indo_dat_treated$treated == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0, ]$education)
)

summary(lm(education ~ treated, data = indo_dat_treated))

# (C-d)
cd <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 0, ]$education)
)
print(cd)

# (C-e)
ce <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 0, ]$education)
)
print(ce)

# (C-f)
cf1 <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 1, ]$education)
)
print(cf1)

cf2 <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 0, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 0, ]$education)
)
print(cf2)

# (C-g)
print(ce - cd)
print(cf2 - cf1)

educ_did <- felm(education ~ treated + program_intensity +
            treated:program_intensity, data = indo_dat_treated)
summary(educ_did)


# (C-j)
wage_did <- felm(log_wage ~ treated + program_intensity +
            treated:program_intensity, data = indo_dat_treated)
summary(wage_did)

## Question 2
qu2_dat <- readRDS("data/df.r")
gamma <- 0.3
t1 <- 5
t2 <- 7
tau1 <- -0.1

# (a)
treatments_by_month <- qu2_dat %>%
    group_by(month, treatment_month) %>%
    summarise(
        avg_outcome = mean(y)
    )

ggplot(
    data = treatments_by_month,
    aes(x = month, y = avg_outcome, group = treatment_month)
) + geom_line(aes(color = treatment_month)) +
    geom_vline(xintercept = t1, color = "red", linetype = "dashed") +
    geom_vline(xintercept = t2, color = "red", linetype = "dashed")

# (c)

# TWFE
qu2_dat$month_from_treat <- ifelse(
    qu2_dat$treatment_month == 0,
    99999,
    qu2_dat$month - qu2_dat$treatment_month
)
qu2_dat$treat <- ifelse(
    qu2_dat$treatment_month == 0,
    0,
    1
)
qu2_dat$month_from_treat <- relevel(
    as.factor(qu2_dat$month_from_treat), ref = "0"
)

twfe <- felm(
    y ~ treat + treat:month_from_treat | month + id,
    data = qu2_dat
)
summary(twfe)

did_results <- data.frame(
    from_treat = -t1:t2,
    twfe_effect = twfe$coefficients[2:14, 1]
)

# Callaway-Sant'anna’s
model <- att_gt(yname = "y", tname = "month", idname = "id",
gname = "treatment_month", data = qu2_dat, panel = T)
w_did_model <- aggte(model, type = "dynamic")

did_treat_5 <- att_gt(yname = "y", tname = "month", idname = "id",
gname = "treatment_month", data = qu2_dat[qu2_dat$treatment_month != t2, ], panel = T)
agg_did_treat_5 <- aggte(did_treat_5, type = "dynamic")

did_treat_7 <- att_gt(yname = "y", tname = "month", idname = "id",
gname = "treatment_month", data = qu2_dat[qu2_dat$treatment_month != t1, ], panel = T)
agg_did_treat_7 <- aggte(did_treat_7, type = "dynamic")

did_results$att_effect <- w_did_model$att.egt
did_results$att_effect_5 <- c(0, 0, agg_did_treat_5$att.egt)
did_results$att_effect_7 <- c(agg_did_treat_7$att.egt, 0, 0)

ggplot(
    data = did_results,
    aes(x = from_treat)
) + geom_line(aes(y = twfe_effect, color = "red", )) +
    geom_line(aes(y = att_effect, color = "blue")) +
    scale_color_manual(
        name="Method", labels = c("TWFE", "CS"), values = c("blue", "red")
    )

# Callaway-Sant'anna’s unbiasness and consistency graph

# real effects
did_results$dgp_effect_early <- did_results$from_treat +
                                tau1 * did_results$from_treat^2
did_results$dgp_effect_later <- gamma * (did_results$from_treat +
                                tau1 * did_results$from_treat^2)
ggplot(
    data = did_results,
    aes(x = from_treat)
) +
    geom_line(aes(y = dgp_effect_early, color = "red")) +
    geom_line(aes(y = dgp_effect_later, color = "blue")) +
    geom_line(aes(y = att_effect_5, color = "green")) +
    geom_line(aes(y = att_effect_7, color = "pink")) +
    scale_color_manual(name = "Legend", labels = c(
        "DGP 5th month effect",
        "DGP 7th month effect",
        "Callaway-Sant'anna’s for 5th month",
        "Callaway-Sant'anna’s for 7th month"
    ), values = c(
        "red", "blue", "green", "pink"
    ))
