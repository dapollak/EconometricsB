library(haven)
library(dplyr)
library(tidyr)
library(lfe)
library(car)
library(stargazer)
library(strucchange)
library(lmtest)
library(pwrss)
library(xtable)

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
for (i in seq(-11, 6)) {
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

ev_study_mmr <- felm(lnm_rate ~ beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
          beta_1 + beta_2 + beta_3 + beta_4 + beta_5 + beta_6 | factor(year) + factor(state):year + disease, 
          data = filter(subset_state_data_1925_1943_long, 1925 <= year & (disease %in% c("mmr", "tb_rate"))))
summary(ev_study_mmr)
stargazer(coef(summary(ev_study_mmr))[, 1:2], flip = TRUE)

ev_study_inf_pne <- felm(lnm_rate ~ beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_1 + beta_2 + beta_3 + beta_4 + beta_5 + beta_6 | factor(year) + factor(state):year + disease,
               data = filter(subset_state_data_1925_1943_long,1925 <= year & (disease %in% c("infl_pneumonia_rate", "tb_rate"))))
summary(ev_study_inf_pne)
stargazer(coef(summary(ev_study_inf_pne))[, 1:2], flip = TRUE)

subset_state_data_1925_1943_long[is.na(subset_state_data_1925_1943_long) | subset_state_data_1925_1943_long == Inf | subset_state_data_1925_1943_long == -Inf] <- NA

ev_study_scarlet <- felm(lnm_rate ~ beta_min_5 + beta_min_4 + beta_min_3 + beta_min_2 + beta_min_1 + beta_0 +
                        beta_1 + beta_2 + beta_3 + beta_4 + beta_5 + beta_6 | factor(year) + factor(state):year + disease,
               data = na.omit(filter(subset_state_data_1925_1943_long, (disease %in% c("scarfever_rate", "tb_rate")))))
summary(ev_study_scarlet)
stargazer(coef(summary(ev_study_scarlet))[, 1:2], flip = TRUE)

stargazer(
  ev_study_mmr,
  ev_study_inf_pne,
  ev_study_scarlet,
  type = "latex"
)

#####################################

# State data fixed effects
state_fe <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37| 0 | disease_year,
             data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))

states_data <- getfe(state_fe)
states_data$diff <- dplyr::lead(c$effect, 48, 1) - c$effect
states_data <- states_data[1:48,]
states_data$state <- (subset_state_data_1925_1943_long %>%
                        select(state) %>%
                        distinct())$state
states_data <- states_data[c("state", "diff")]

# state data seperate regression (Table A.1 & Table A.4 & Figure 9)
e_all <- c()
e_significant <- c()
needed_sample <- c()
se_vals <- c()
p_vals <- c()
for (state_s in (subset_state_data_1925_1943_long %>% select(state) %>% distinct())$state) {
  filtered_data <- filter(
    subset_state_data_1925_1943_long,
    state == state_s & disease %in% c("mmr", "tb_rate")
  )

  reg <- lm(lnm_rate ~ treated:post37 + treated:year_c + treated +
  year_c + post37, data = filtered_data)

  lh <- linearHypothesis(reg, c(
    "treatedTRUE:post37TRUE=0"
  ), white.adjust = "hc1")

  effect <- reg$coefficients["treatedTRUE:post37TRUE"]

  # power analysis
  p <- pwrss.t.reg(
    beta1 = effect,
    # sdy = sd(filtered_data$lnm_rate, na.rm = T),
    sdx = sqrt(0.5 * (1 - 0.5)), #summary(reg)$coefficients[5, 2]^2,
    k = 5,
    r2 = 0.5,
    power = .7,
    alpha = 0.05,
    alternative = "not equal"
  )

  needed_sample <- c(needed_sample, p$n)

  # store effect and whether significant
  e_all <- c(e_all, effect)
  if (lh$"Pr(>F)"[2] < 0.05) {
    e_significant <- c(e_significant, effect)
  } else {
    e_significant <- c(e_significant, 0)
  }

  se_vals <- c(se_vals, summary(reg)$coefficients["treatedTRUE:post37TRUE","Std. Error"])
  p_vals <- c(p_vals, summary(reg)$coefficients["treatedTRUE:post37TRUE","Pr(>|t|)"])
}
states_data$seperate_effects_all <- e_all
states_data$seperate_effects_significant <- e_significant
states_data$power_sample_size <- needed_sample
states_data$se <- se_vals
states_data$p_vals <- p_vals
rownames(states_data) <- states_data$state

# States Effect (Table A.1, has data for Figure 3)
states_effects <- states_data[c(
  "seperate_effects_all",
  "seperate_effects_significant",
  "se",
  "p_vals"
)]

# dump csv for the heatmap
write.csv(states_data, "/tmp/heatmap_mmr.csv")

# Table A.4
power_table <- states_data[c("seperate_effects_all", "power_sample_size")]
power_table <- power_table[order(power_table$power_sample_size), ]
writeLines(xtable(power_table), "power.tex")


## Table 3 & Figure 2 - National trend year break paper replication
national_data$all_break_output <- log(national_data$all_tot) -
                              dplyr::lag(log(national_data$all_tot))
national_data$mmr_break_output <- log(national_data$mmr) -
                              dplyr::lag(log(national_data$mmr))
national_data$tuberculosis_break_output <- log(national_data$tuberculosis_total) -
                              dplyr::lag(log(national_data$tuberculosis_total))
national_data$inf_pne_break_output <- log(national_data$influenza_pneumonia_total) -
                              dplyr::lag(log(national_data$influenza_pneumonia_total))
national_data$scarlet_break_output <- log(national_data$scarlet_fever_tot) -
                              dplyr::lag(log(national_data$scarlet_fever_tot))

print("all,mmr,inf_pne,scarlet,tuberculosis")
nation_trend_break <- data.frame()
for (tau in 1933:1942) {
  all_break_reg <- lm(all_break_output ~ I(year >= tau), data = national_data)
  all_t_val <- coeftest(all_break_reg, vcov.=NeweyWest(all_break_reg, lag=1, adjust=T))[6]
  
  mmr_break_reg <- lm(mmr_break_output ~ I(year >= tau), data = national_data)
  mmr_t_val <- coeftest(mmr_break_reg, vcov.=NeweyWest(mmr_break_reg, lag=1, adjust=T))[6]

  inf_pne_break_reg <- lm(inf_pne_break_output ~ I(year >= tau), data = national_data)
  inf_pne_t_val <- coeftest(inf_pne_break_reg, vcov.=NeweyWest(inf_pne_break_reg, lag=1, adjust=T))[6]

  scarlet_break_reg <- lm(scarlet_break_output ~ I(year >= tau), data = national_data)
  scarlet_t_val <- coeftest(scarlet_break_reg, vcov.=NeweyWest(scarlet_break_reg, lag=1, adjust=T))[6]

  tub_break_reg <- lm(tuberculosis_break_output ~ I(year >= tau), data = national_data)
  tub_t_val <- coeftest(tub_break_reg, vcov.=NeweyWest(tub_break_reg, lag=1, adjust=T))[6]

  print(sprintf("%s,%s,%s,%s,%s,%s", tau, all_t_val^2, mmr_t_val^2, inf_pne_t_val^2, scarlet_t_val^2, tub_t_val^2))
  nation_trend_break <- rbind(nation_trend_break, data.frame(
    year=tau, all_f_val=all_t_val^2, mmr_f_val=mmr_t_val^2, inf_pne_f_val=inf_pne_t_val^2,
    scarlet_f_val=scarlet_t_val^2, tub_f_val=tub_t_val^2
  ))
}

# Table 3
replication_year <- data.frame(
  year=c(
    nation_trend_break[which.max(nation_trend_break$all_f_val),]$year,
    nation_trend_break[which.max(nation_trend_break$mmr_f_val),]$year,
    nation_trend_break[which.max(nation_trend_break$inf_pne_f_val),]$year,
    nation_trend_break[which.max(nation_trend_break$scarlet_f_val),]$year,
    nation_trend_break[which.max(nation_trend_break$tub_f_val),]$year
  ),
  f_stat=c(
    nation_trend_break[which.max(nation_trend_break$all_f_val),]$all_f_val,
    nation_trend_break[which.max(nation_trend_break$mmr_f_val),]$mmr_f_val,
    nation_trend_break[which.max(nation_trend_break$inf_pne_f_val),]$inf_pne_f_val,
    nation_trend_break[which.max(nation_trend_break$scarlet_f_val),]$scarlet_f_val,
    nation_trend_break[which.max(nation_trend_break$tub_f_val),]$tub_f_val
  )
)
rownames(replication_year) <- c("All", "MMR", "Influenza", "Scarlet", "Tuberculosis")
writeLines(xtable(states_trend_break_max), "trend_break_table.tex")

## Table A.2 - States trend year break extension
subset_state_data_1925_1943_long$break_output <- subset_state_data_1925_1943_long$lnm_rate - dplyr::lag(subset_state_data_1925_1943_long$lnm_rate)

states_trend_break_full <- data.frame()
for (state_s in subset_state_data_1925_1943_long$state[!duplicated(subset_state_data_1925_1943_long$state)]) {
  for (tau in 1933:1940) {
    mmr_break_reg <- lm(
      break_output ~ I(year >= tau),
      data = filter(subset_state_data_1925_1943_long, state == state_s & disease == "mmr")
    )
    tryCatch(
      {
        mmr_t_val <- coeftest(mmr_break_reg, vcov.=NeweyWest(mmr_break_reg, lag=1, adjust=T))[6]

        states_trend_break_full <- rbind(
          states_trend_break_full,
          data.frame(state=state_s, year=tau, f_stat=mmr_t_val^2)
        )
      },
      error = function(e) { }
    )
  }
}

states_trend_break_max <- states_trend_break_full %>%
            group_by(state) %>%
            mutate(f_stat_max=max(f_stat)) %>%
            filter(f_stat == f_stat_max)

# Table A.2
states_trend_break_max <- subset(states_trend_break_max, select = -c(f_stat_max))
writeLines(xtable(states_trend_break_max), "states_year_break.tex")