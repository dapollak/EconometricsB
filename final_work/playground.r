library(haven)
library(dplyr)
library(tidyr)

rm(list = ls())
options(scipen = 999)

sf_mortality_data <- read_dta("final_work/data/datasets/20080229_SFmortality.dta")
sf_mortality_race_data <- read_dta("final_work/data/datasets/20080229_SFmortality_Race.dta")
state_data <- read_dta("final_work/data/datasets/20080229_state_data.dta")
state_race_data <- read_dta("final_work/data/datasets/20080229_state_race_data.dta")
national_data <- read_dta("final_work/data/datasets/20080229_national_data.dta")

### Event study state level
# Filter the data to keep only the years between 1925 and 1943
state_data_1925_1943 <- state_data %>%
  filter(year >= 1925 & year <= 1943)

subset_state_data_1925_1943 <- state_data_1925_1943[, c("state", "year", "mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate")]

#Reshape long
subset_state_data_1925_1943_long <- subset_state_data_1925_1943 %>%
  pivot_longer(cols = c("mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate"),
               names_to = "disease",
               values_to = "m_rate") 

# Create new variables
subset_state_data_1925_1943_long <- subset_state_data_1925_1943_long %>%
  mutate(treated = (disease == "mmr" | disease == "infl_pneumonia_rate" | disease == "scarfever_rate"),
         post37 = (year >= 1937),
         year_c = year - 1937,
         lnm_rate = log(m_rate))


subset_state_data_1925_1943_long$state_post37 <- as.numeric(interaction(subset_state_data_1925_1943_long$state, subset_state_data_1925_1943_long$post37))  
subset_state_data_1925_1943_long$disease_year <- as.numeric(interaction(subset_state_data_1925_1943_long$disease, subset_state_data_1925_1943_long$year))  

# add lags and leads
for (i in seq(-5, 6)) {
    if (i < 0) {
        col_name <- sprintf("beta_min_%s", abs(i))
    }
    else {
        col_name <- sprintf("beta_%s", i)
    }
    subset_state_data_1925_1943_long[col_name] <-
                            ifelse(subset_state_data_1925_1943_long$year_c - i == 0 & 
                                    subset_state_data_1925_1943_long$treated, 1, 0)
}

ev_study1 <- lm(lnm_rate ~ treated + beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_6 + beta_5 + beta_4 + beta_3 + beta_2 + beta_1, 
               data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
summary(ev_study1)

ev_study2 <- lm(lnm_rate ~ treated + beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_6 + beta_5 + beta_4 + beta_3 + beta_2 + beta_1,
               data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
summary(ev_study2)

ev_study3 <- lm(lnm_rate ~ treated + beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_6 + beta_5 + beta_4 + beta_3 + beta_2 + beta_1,
               data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
summary(ev_study3)
