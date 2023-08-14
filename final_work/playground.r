library(haven)
library(dplyr)
library(tidyr)
library(lfe)
library(car)

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
                                    (0 | subset_state_data_1925_1943_long$treated), 1, 0)
}

ev_study1_f <- felm(lnm_rate ~ beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_6 + beta_5 + beta_4 + beta_3 + beta_2 + beta_1 | factor(state):year + disease, 
               data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "infl_pneumonia_rate", "tb_rate")))
summary(ev_study1_f)

linearHypothesis(ev_study1_f, c(
  "beta_min_5=0",
  "beta_min_4=0",
  "beta_min_3=0",
  "beta_min_2=0",
  "beta_min_1=0"
), white.adjust = "hc1")

ev_study2_f <- felm(lnm_rate ~ beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_6 + beta_5 + beta_4 + beta_3 + beta_2 + beta_1 | factor(state):year + disease,
               data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
summary(ev_study2_f)

linearHypothesis(ev_study2_f, c(
  "beta_min_5=0",
  "beta_min_4=0",
  "beta_min_3=0",
  "beta_min_2=0",
  "beta_min_1=0"
), white.adjust = "hc1")

# State data fixed effects
state_fe <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37| 0 | disease_year,
             data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))

fixed_effects <- getfe(state_fe)
fixed_effects$diff <- lead(c$effect, 48, 1) - c$effect
fixed_effects <- fixed_effects[1:48,]
fixed_effects$state <- (subset_state_data_1925_1943_long %>%
                        select(state) %>%
                        distinct())$state
fixed_effects <- fixed_effects[c("state", "diff")]

# state data seperate regression
e_all <- c()
e_significant <- c()
for (state_s in (subset_state_data_1925_1943_long %>% select(state) %>% distinct())$state) {
  reg <- lm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_state_data_1925_1943_long, state == state_s & disease %in% c("mmr", "tb_rate")))
  
  lh <- linearHypothesis(reg, c(
    "treatedTRUE:post37TRUE=0"
  ), white.adjust = "hc1")

  e_all <- c(e_all, reg$coefficients["treatedTRUE:post37TRUE"])
  if (lh$"Pr(>F)"[2] < 0.05) {
    e_significant <- c(e_significant, reg$coefficients["treatedTRUE:post37TRUE"])
  } else {
    e_significant <- c(e_significant, 0)
  }
}
fixed_effects$seperate_effects_all <- e_all
fixed_effects$seperate_effects_significant <- e_significant

write.csv(fixed_effects, "/tmp/hm.csv")
