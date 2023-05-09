library(dplyr)
library(plm)
library(lfe)
library(haven)
options(scipen = 999)

indo_dat <- read_dta("data/Indo_Schooling.dta")

## (B)

reg1 <- lm(log_wage ~ education, data = indo_dat)

## (C)
intense_levels <- indo_dat[c("program_intensity", "num_schools", "education")]

# (C-a)
intense_levels <- intense_levels %>%
    group_by(program_intensity) %>%
    summarise(
        num_schools_avg = mean(num_schools),
        education_avg = mean(education)
    )

# (C-b)
print(
    mean(indo_dat[indo_dat$program_intensity == 1, ]$education) -
    mean(indo_dat[indo_dat$program_intensity == 0, ]$education)
)

# (C-c)
indo_dat_treated <- indo_dat %>%
    mutate(
        treated = ifelse(birth_year < 62, 0, 1)
    )

print(
    mean(indo_dat_treated[indo_dat_treated$treated == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0, ]$education)
)

# (C-d)
cd <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 0, ]$education)
)

# (C-e)
ce <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 0, ]$education)
)

# (C-f)
cf1 <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 1, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 1, ]$education)
)
cf2 <- (
    mean(indo_dat_treated[indo_dat_treated$treated == 1 &
        indo_dat_treated$program_intensity == 0, ]$education) -
    mean(indo_dat_treated[indo_dat_treated$treated == 0 &
        indo_dat_treated$program_intensity == 0, ]$education)
)

did1 <- felm(education ~ treated + program_intensity + treated:program_intensity, data = indo_dat_treated)


# (C-j)
wage_did <- felm(log_wage ~ treated + program_intensity + treated:program_intensity, data = indo_dat_treated)


## Question 2