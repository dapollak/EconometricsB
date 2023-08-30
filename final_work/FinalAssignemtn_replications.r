
#########   - Section 1- Replication   ##########
# The number of the tables and figures corresponds to the numbers in the assignment and not the original paper 

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



setwd("G:/האחסון שלי/לימודים/אקונומטריקה/אקונומטריקה ב/פרויקט גמר/AEJApp_Datasets")

### Sub section- Descriptive Statistics

## Table 1- National Mortality Statistics (deaths per 100,000), Mean

national_data <- read_dta("national_data.dta")

national_data_1920 <- subset(national_data, year == 1920)
national_data_1950 <- subset(national_data, year == 1950)
national_data_1925_1936 <- subset(national_data, year >= 1925 & year <= 1936)
national_data_1937_1943 <- subset(national_data, year >= 1937 & year <= 1943)


national_summary_1920 <- national_data_1920 %>%
  summarize_at(vars(all_tot, mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total, diabetes_total, heart_total, cancer_total), mean)

national_summary_1950 <- national_data_1950 %>%
  summarize_at(vars(all_tot, mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total, diabetes_total, heart_total, cancer_total), mean)

national_summary_1925_1936 <- national_data_1925_1936 %>%
  summarize_at(vars(all_tot, mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total, diabetes_total, heart_total, cancer_total), mean)

national_summary_1937_1943 <- national_data_1937_1943 %>%
  summarize_at(vars(all_tot, mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total, diabetes_total, heart_total, cancer_total), mean)


all_summaries <- rbind(national_summary_1920, national_summary_1950, national_summary_1925_1936, national_summary_1937_1943)
transposed_table <- t(all_summaries)
column_names <- c("1920", "1950", "1925-1936", "1937-1943")
colnames(transposed_table) <- column_names
latex_table <- xtable(transposed_table, digits = 2)
file_path <- "replication_table1_panel_a.tex"
print(latex_table, file = file_path, floating = FALSE)



## Table 2- State Level Mortality Statistics (deaths per 100,000), Mean

state_data <- read_dta("state_data.dta")

state_data_1925_1936 <- subset(state_data, year >= 1925 & year <= 1936)
state_data_1937_1943 <- subset(state_data, year >= 1937 & year <= 1943)

state_summary_1925_1936 <- state_data_1925_1936 %>% 
  summarize_at(vars(mmr, infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), ~ mean(., na.rm = TRUE)) %>%
  mutate_at(vars(infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), ~ . * 100000)
state_summary_1937_1943 <- state_data_1937_1943 %>% 
  summarize_at(vars(mmr, infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), ~ mean(., na.rm = TRUE)) %>%
  mutate_at(vars(infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), ~ . * 100000)

all_summaries <- rbind(state_summary_1925_1936, state_summary_1937_1943)
transposed_table <- t(all_summaries)
column_names <- c("1925-1936", "1937-1943")
colnames(transposed_table) <- column_names
latex_table <- xtable(transposed_table)
file_path <- "replication_table1_panel_b.tex"
print(latex_table, file = file_path, floating = FALSE)



### Sub section-  Graphical and Formal test for Break Year in Mortality Trends

## Figure 1- Mortality Rates

national_data <- read_dta("national_data.dta")
subset_national_data <- national_data[, c("mmr", "influenza_pneumonia_total", "scarlet_fever_tot", "tuberculosis_total", "year")]

subset_national_data_long <- subset_national_data %>%
  pivot_longer(cols = c(mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total),
               names_to = "disease",
               values_to = "m_rate") 

subset_national_data_long <- subset_national_data_long %>%
  mutate(lnm_rate = log(m_rate))


plot1<-ggplot(data = subset_national_data_long, aes(x = year, y = lnm_rate, color = disease)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 1937, linetype = "dashed", color = "red") +
  labs(x = "Year",
       y = "Log Mortality Rate",
       color = "Disease")+
  scale_color_manual(
    values = c("black", "green", "blue", "red"),
    labels = c("Influenza and Pneumonia", "Maternal Mortality", "Scarlet Fever", "Tuberculosis")
  )+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        legend.margin =margin(10,10,10,10),
        axis.line = element_line(color = "black"))+
  scale_x_continuous(breaks = seq(1920, 1950, by = 5), limits = c(1920, 1950)) +
  scale_y_continuous(breaks = seq(-3, 7.5, by = 1.5), limits = c(-3, 7)) +
  geom_hline(yintercept = seq(-3, 7.5, by = 1.5), linetype = "dashed", color = "lightgray")

print(plot1)

ggsave("Log Mortality Rate.png", plot1, width = 8, height = 4)



## Table 3- Trend break year test results
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
writeLines(print(xtable(replication_year)), "trend_break_table.tex")

## Figure 2- Trend break year test results
write.csv(nation_trend_break, "break_year.csv")

### Sub section- Main results

## Table 4 - Replication of table 4 Panel A

national_data <- read_dta("national_data.dta")
national_data_1925_1943 <- national_data %>%
  filter(year >= 1925 & year <= 1943)

subset_national_data_1925_1943 <- national_data_1925_1943[, c("mmr", "influenza_pneumonia_total", "scarlet_fever_tot", "tuberculosis_total", "year")]

subset_national_data_1925_1943_long <- subset_national_data_1925_1943 %>%
  pivot_longer(cols = c(mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total),
               names_to = "disease",
               values_to = "m_rate") 

subset_national_data_1925_1943_long <- subset_national_data_1925_1943_long %>%
  mutate(treated = (disease == "mmr" | disease == "influenza_pneumonia_total" | disease == "scarlet_fever_tot"),
         post37 = (year >= 1937),
         year_c = year - 1937,
         lnm_rate = log(m_rate))


reg1 <-feols(lnm_rate~ treated:post37 + treated:year_c + treated + year_c + post37, se="hetero",
             data = filter(subset_national_data_1925_1943_long, disease %in% c("mmr", "tuberculosis_total")))
reg2 <-feols(lnm_rate~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37, se="hetero",
             data = filter(subset_national_data_1925_1943_long, disease %in% c("mmr", "tuberculosis_total")))

reg3 <-feols(lnm_rate~ treated:post37 + treated:year_c + treated + year_c + post37, se="hetero",
             data = filter(subset_national_data_1925_1943_long, disease %in% c("influenza_pneumonia_total", "tuberculosis_total")))
reg4 <-feols(lnm_rate~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37, se="hetero",
             data = filter(subset_national_data_1925_1943_long, disease %in% c("influenza_pneumonia_total", "tuberculosis_total")))

reg5 <-feols(lnm_rate~ treated:post37 + treated:year_c + treated + year_c + post37, se="hetero",
             data = filter(subset_national_data_1925_1943_long, disease %in% c("scarlet_fever_tot", "tuberculosis_total")))
reg6 <-feols(lnm_rate~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37, se="hetero",
             data = filter(subset_national_data_1925_1943_long, disease %in% c("scarlet_fever_tot", "tuberculosis_total")))

stargazer(reg1,reg2, reg3, reg4, reg5, reg6, type = "latex")
table4_panel_a <- list(reg1, reg2, reg3, reg4, reg5, reg6)
latex_code <- texreg(table4_panel_a)
writeLines(latex_code, "table4_panel_a.tex")


## Table 5 - Replication of table 4 panel B

state_data <- read_dta("state_data.dta")
state_data_1925_1943 <- state_data %>%
  filter(year >= 1925 & year <= 1943)

subset_state_data_1925_1943 <- state_data_1925_1943[, c("state", "year", "mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate")]

subset_state_data_1925_1943_long <- subset_state_data_1925_1943 %>%
  pivot_longer(cols = c("mmr", "infl_pneumonia_rate", "scarfever_rate", "tb_rate"),
               names_to = "disease",
               values_to = "m_rate")

subset_state_data_1925_1943_long <- subset_state_data_1925_1943_long %>%
  mutate(treated = (disease == "mmr" | disease == "infl_pneumonia_rate" | disease == "scarfever_rate"),
         post37 = (year >= 1937),
         year_c = year - 1937,
         lnm_rate = log(m_rate))

subset_state_data_1925_1943_long$disease_year <- as.numeric(interaction(subset_state_data_1925_1943_long$disease, subset_state_data_1925_1943_long$year))  
subset_state_data_1925_1943_long$state_post37 <- as.numeric(interaction(subset_state_data_1925_1943_long$state, subset_state_data_1925_1943_long$post37))  


reg7 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c |state_post37|0|disease_year, 
             data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
reg8 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + state_post37:year_c |state_post37|0|disease_year, 
              data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))

reg9 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c |state_post37|0|disease_year, 
           data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
reg10 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + state_post37:year_c |state_post37|0|disease_year, 
            data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))

subset_state_data_1925_1943_long <- subset(subset_state_data_1925_1943_long, m_rate != 0)
reg11 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c |state_post37|0|disease_year, 
            data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
reg12 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + state_post37:year_c |state_post37|0|disease_year, 
            data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))

table4_panel_b <- list(reg7, reg8, reg9, reg10, reg11, reg12)
latex_code <- texreg(table4_panel_b)
writeLines(latex_code, "table4_panel_b.tex")
