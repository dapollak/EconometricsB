
#########   - Section 3- Extensions   ##########
# The number of the tables and figures corresponds to the numbers in the assignment and not the original paper 
# After producing the relevant data the heat maps produced using the website: http://www.heatmapper.ca/geomap/
# Few final tables were edited in overleaf, but all the coefficients were produced using r. 


## Clear desk
rm(list=ls())

## Garbage collection
gc()

## Avoid scientific notation
options(scipen=999)

## Load libraries  
library(dplyr)
library(plm)
library(lfe)
library(haven)
library(ggplot2)
library(did)
library(coefplot)
library(stargazer)
library(tidyverse)
library(texreg)
library(xtable)
library(strucchange)
library(lmtest)
library(fixest)
library(pwr)
library(readxl)
library(rgdal)
library(sf)
library(leaflet)
library(writexl)



setwd("G:/האחסון שלי/לימודים/אקונומטריקה/אקונומטריקה ב/פרויקט גמר/AEJApp_Datasets")

### Sub section- State Heterogeneity

## Figure 3 and Table A.1 in the appendix- Treatment effects by state
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

# States Effect (Table A.1)
states_effects <- states_data[c(
  "seperate_effects_all",
  "seperate_effects_significant",
  "se",
  "p_vals"
)]
writeLines(xtable(states_effects), "table_did_weighted.tex")

# dump csv for the heatmap (figure 3)
write.csv(states_data, "heatmap_mmr.csv")

## Table 6- Extreme State Effects
extreme_states_effects <- states_effects[
  c(
    which.max(states_effects$seperate_effects_all),
    which.min(states_effects$seperate_effects_all)
  ),
][c("seperate_effects_all")]
print(xtable(extreme_states_effects))


## Figure 4 - Black Population Ratio
#The data is from census 1930 and attached in excel file. Heat maps produced using the website: http://www.heatmapper.ca/geomap/ 


## Figure 5
#from website: https://en.wikipedia.org/wiki/75th_United_States_Congress


### Sub section- Diff-in-Diff Model 

## Table A.2 in the appendix- Trend break year by states
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

states_trend_break_max <- subset(states_trend_break_max, select = -c(f_stat_max))
writeLines(xtable(states_trend_break_max), "states_year_break.tex")

## Figure 6 -8 + Table A.3 - Event study

state_data <- read_dta("state_data.dta")
state_data_1925_1943 <- state_data %>%
  filter(year >= 1925 & year <= 1943)

subset_state_data_1925_1943 <- state_data_1925_1943[, c("state", "year", "mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate")]

subset_state_data_1925_1943_long <- subset_state_data_1925_1943 %>%
  pivot_longer(cols = c("mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate"),
               names_to = "disease",
               values_to = "m_rate") 

subset_state_data_1925_1943_long <- subset_state_data_1925_1943_long %>%
  mutate(treated = as.numeric(disease == "mmr" | disease == "infl_pneumonia_rate" | disease == "scarfever_rate"),
         year_c = year - 1937,
         lnm_rate = log(m_rate))

ev_study_mmr <- feols(lnm_rate~ i(year,treated,ref=1937) | year + treated + state^year, se="hetero",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
plt1 = coefplot(ev_study_mmr, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.5,1.5))
summary(ev_study_mmr)

ev_study_inf_pne <- feols(lnm_rate~ i(year,treated,ref=1937) | year + treated + state^year, se="hetero",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
plt2 = coefplot(ev_study_inf_pne, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.5,1.5))
summary(ev_study_inf_pne)

ev_study_scarlet <- feols(lnm_rate~ i(year,treated,ref=1937) | year + treated + state^year, se="hetero",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
plt3 = coefplot(ev_study_scarlet, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.5,1.5))
summary(ev_study_scarlet)


## Table 7- Weighted Diff-in-Diff estimation
census_30 <- "census_1930.xlsx"
census_30 <- read_excel(census_30)

state_data <- read_dta("state_data.dta")

state_data_1925_1943 <- state_data %>%
  filter(year >= 1925 & year <= 1943)

subset_state_data_1925_1943 <- state_data_1925_1943[, c("state", "year", "mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate")]

subset_state_data_1925_1943_long <- subset_state_data_1925_1943 %>%
  pivot_longer(cols = c("mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate"),
               names_to = "disease",
               values_to = "m_rate") 

merged_state_data <- left_join(subset_state_data_1925_1943_long, census_30, by = "state")


merged_state_data <- merged_state_data %>%
  mutate(treated = (disease == "mmr" | disease == "infl_pneumonia_rate" | disease == "scarfever_rate"),
         post37 = (year >= 1937),
         year_c = year - 1937,
         lnm_rate = log(m_rate)
         )

merged_state_data$disease_year <- as.numeric(interaction(merged_state_data$disease, merged_state_data$year))  
merged_state_data$state_post37 <- as.numeric(interaction(merged_state_data$state, merged_state_data$post37))  
merged_state_data$Population_total <- as.numeric(merged_state_data$Population_total)  

merged_state_data <- subset(merged_state_data, m_rate != 0)

mmr <- merged_state_data %>%  filter(disease %in% c("mmr", "tb_rate"))
reg13 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c |state_post37|0|disease_year, 
             data = mmr,
             weights=mmr$Population_total)

reg14 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + state_post37:year_c |state_post37|0|disease_year, 
              data = mmr,
              weights=mmr$Population_total)


infl_pne <- merged_state_data %>%  filter(disease %in% c("infl_pneumonia_rate", "tb_rate"))
reg15 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c |state_post37|0|disease_year, 
             data = infl_pne,
             weights=infl_pne$Population_total)

reg16 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + state_post37:year_c |state_post37|0|disease_year, 
              data = infl_pne,
              weights=infl_pne$Population_total)


merged_state_data <- subset(merged_state_data, m_rate != 0)
scarfever <- merged_state_data %>%  filter(disease %in% c("scarfever_rate", "tb_rate"))

reg17 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c |state_post37|0|disease_year, 
              data = scarfever,
              weights=scarfever$Population_total)

reg18 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + state_post37:year_c |state_post37|0|disease_year, 
              data = scarfever,
              weights=scarfever$Population_total)


table_did_weighted <- list(reg13, reg14, reg15, reg16, reg17, reg18)
latex_code <- texreg(table_did_weighted)
writeLines(latex_code, "table_did_weighted.tex")



### Sub section-Power Calculations

## Table A.4 in the appendix- Power Calculations
power_table <- states_data[c("seperate_effects_all", "power_sample_size")]
power_table <- power_table[order(power_table$power_sample_size), ]
writeLines(xtable(power_table), "power.tex")