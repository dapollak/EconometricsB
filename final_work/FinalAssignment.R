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



# setwd("/Users/danielpollak/Library/CloudStorage/GoogleDrive-pollak.daniel@gmail.com/My Drive/University/Economics/Econometrics B/final_work")


###Replication of table 2- structural break test###

national_data <- read_dta("final_work/data/datasets/20080229_national_data.dta")
national_data$lnmmr <- log(national_data$mmr)
national_data$lnpne <- log(national_data$influenza_pneumonia_total)
national_data$lnsc <- log(national_data$scarlet_fever_tot)
national_data$lntb <- log(national_data$tuberculosis_total)
national_data$lndia <- log(national_data$diabetes_total)
national_data$lnheart <- log(national_data$heart_total)
national_data$lnall <- log(national_data$all_tot)
national_data$lncancer <- log(national_data$cancer_total)


national_data_ts <- ts(national_data$lnmmr, start = national_data$year[1])



# Set up variables to store the results
sig <- rep(NA, 10)
fstat <- rep(NA, 10)
maxf <- NA

# Perform the loop to test for breakpoints
for (i in 33:42) {
  y_var <- paste0("y19", i)
  national_data[[y_var]] <- ifelse(national_data$year >= 1900 + i, 1, 0)
  fit <- lm(national_data$lnmmr ~ national_data[[y_var]])
  national_data$fstat[national_data$year == (1900 + i)] <- summary(fit)$fstatistic[1]
  test_result <- coef(summary(fit))["national_data[[y_var]]", "Pr(>|t|)"]
  national_data$sig[national_data$year == (1900 + i)] <- test_result
}

# Identify the maximum fstat value
maxf <- max(national_data$fstat, na.rm = TRUE)

# Identify the breakpoint based on the maximum fstat value
break_year <- national_data$year[national_data$fstat == maxf]

# Print the results
print(national_data[c("lnmmr", "year", "break", "fstat", "sig")])






### Replication of table 1- summary statistics
##Panel A
national_data <- read_dta("final_work/data/datasets/20080229_national_data.dta")

#average national mortality rate- all
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
latex_table <- xtable(transposed_table)
file_path <- "replication_table1_panel_a.tex"
print(latex_table, file = file_path, floating = FALSE)


#average national mortality rate- by race
national_summary_1920_r <- national_data_1920 %>% 
  summarize_at(vars(all_w, all_nw, mmr_w, mmr_nw, influenza_pneumonia_w, influenza_pneumonia_nw, scarlet_fever_w, scarlet_fever_nw, tuberculosis_w, tuberculosis_nw), mean)
national_summary_1950_r <- national_data_1950 %>% 
  summarize_at(vars(all_w, all_nw, mmr_w, mmr_nw, influenza_pneumonia_w, influenza_pneumonia_nw, scarlet_fever_w, scarlet_fever_nw, tuberculosis_w, tuberculosis_nw), mean)
national_summary_1925_1936_r <- national_data_1925_1936 %>% 
  summarize_at(vars(all_w, all_nw, mmr_w, mmr_nw, influenza_pneumonia_w, influenza_pneumonia_nw, scarlet_fever_w, scarlet_fever_nw, tuberculosis_w, tuberculosis_nw), mean)
national_summary_1937_1943_r <- national_data_1937_1943 %>% 
  summarize_at(vars(all_w, all_nw, mmr_w, mmr_nw, influenza_pneumonia_w, influenza_pneumonia_nw, scarlet_fever_w, scarlet_fever_nw, tuberculosis_w, tuberculosis_nw), mean)


all_summaries <- rbind(national_summary_1920_r, national_summary_1950_r, national_summary_1925_1936_r, national_summary_1937_1943_r)
transposed_table <- t(all_summaries)
column_names <- c("1920", "1950", "1925-1936", "1937-1943")
colnames(transposed_table) <- column_names
latex_table <- xtable(transposed_table)
file_path <- "replication_table1_panel_a_race.tex"
print(latex_table, file = file_path, floating = FALSE)



##Panel B
state_data <- read_dta("data/datasets/state_data.dta")

#average state mortality rate- all

state_data_1925_1936 <- subset(state_data, year >= 1925 & year <= 1936)
state_data_1937_1943 <- subset(state_data, year >= 1937 & year <= 1943)

state_summary_1925_1936 <- state_data_1925_1936 %>% 
  summarize_at(vars(mmr, infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), mean) %>%
  mutate_at(vars(infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), ~ . * 100000)
state_summary_1937_1943 <- state_data_1937_1943 %>% 
  summarize_at(vars(mmr, infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), mean) %>%
  mutate_at(vars(infl_pneumonia_rate, scarfever_rate, tb_rate, diabetes_rate, heartd_rate, cancer_rate), ~ . * 100000)


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


#average state mortality rate- by race
state_data_race <- read_dta("data/datasets/state_data_race.dta")

state_data_race <- state_data_race %>% arrange(race)

summary_1925_1936_state_r <- state_data_race %>%
  filter(year >= 1925 & year <= 1936) %>%
  group_by(race) %>%
  summarize(
    mean_mmr = mean(mmr, na.rm = TRUE),
    mean_flupneu_rate = mean(flupneu_rate, na.rm = TRUE),
    mean_sf_rate = mean(sf_rate, na.rm = TRUE),
    mean_tb_rate = mean(tb_rate, na.rm = TRUE)
  )

summary_1937_1943_state_r <- state_data_race %>%
  filter(year >= 1937 & year <= 1943) %>%
  group_by(race) %>%
  summarize(
    mean_mmr = mean(mmr, na.rm = TRUE),
    mean_flupneu_rate = mean(flupneu_rate, na.rm = TRUE),
    mean_sf_rate = mean(sf_rate, na.rm = TRUE),
    mean_tb_rate = mean(tb_rate, na.rm = TRUE)
  )


all_summaries <- rbind(summary_1925_1936_state_r, summary_1937_1943_state_r)
transposed_table <- t(all_summaries)
column_names <- c("1925-1936 others", "1925-1936 white", "1937-1943 others", "1937-1943 white" )
colnames(transposed_table) <- column_names
latex_table <- xtable(transposed_table)
file_path <- "replication_table1_panel_b_race.tex"
print(latex_table, file = file_path, floating = FALSE)


print(summary_1925_1936_state_r)
print(summary_1937_1943_state_r)


###Replication of table 4- Diff-in-Diff results

##Panel A- National-level data, all years, 1925–1943
national_data <- read_dta("final_work/data/datasets/20080229_national_data.dta")
# Filter the data to keep only the years between 1925 and 1943
national_data_1925_1943 <- national_data %>%
  filter(year >= 1925 & year <= 1943)

subset_national_data_1925_1943 <- national_data_1925_1943[, c("mmr", "influenza_pneumonia_total", "scarlet_fever_tot", "tuberculosis_total", "year")]

#Reshape long
subset_national_data_1925_1943_long <- subset_national_data_1925_1943 %>%
  pivot_longer(cols = c(mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total),
               names_to = "disease",
               values_to = "m_rate") 

# Create new variables
subset_national_data_1925_1943_long <- subset_national_data_1925_1943_long %>%
  mutate(treated = (disease == "mmr" | disease == "influenza_pneumonia_total" | disease == "scarlet_fever_tot"),
         post37 = (year >= 1937),
         year_c = year - 1937,
         lnm_rate = log(m_rate))


# Run regression for each disease

reg1 <- lm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_national_data_1925_1943_long, disease %in% c("mmr", "tuberculosis_total")))

p <- pwrss.t.reg(
  beta1 = reg1$coefficients[5],
  # sdy = sd(filtered_data$lnm_rate, na.rm = T),
  sdx = sqrt(0.5 * (1 - 0.5)), #summary(reg)$coefficients[5, 2]^2,
  k = 5,
  r2 = 0.5,
  power = .7,
  alpha = 0.05,
  alternative = "not equal"
)

reg2 <- lm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_national_data_1925_1943_long, disease %in% c("mmr", "tuberculosis_total")))
reg3 <- lm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_national_data_1925_1943_long, disease %in% c("influenza_pneumonia_total", "tuberculosis_total")))
reg4 <- lm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_national_data_1925_1943_long, disease %in% c("influenza_pneumonia_total", "tuberculosis_total")))
reg5 <- lm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_national_data_1925_1943_long, disease %in% c("scarlet_fever_tot", "tuberculosis_total")))
reg6 <- lm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37, 
               data = filter(subset_national_data_1925_1943_long, disease %in% c("scarlet_fever_tot", "tuberculosis_total")))

stargazer(reg1,reg2, reg3, reg4, reg5, reg6, type = "latex")
table4_panel_a <- list(reg1, reg2, reg3, reg4, reg5, reg6)
latex_code <- texreg(table4_panel_a)
writeLines(latex_code, "table4_panel_a.tex")

##Panel B- State-level, all years, 1925–1943
state_data <- read_dta("final_work/data/datasets/20080229_state_data.dta")
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


# Run regression for each disease

reg7 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37|0|disease_year, 
           data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
reg8 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37 + state_post37:year_c |state_post37|0|disease_year, 
           data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
reg9 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37|0|disease_year, 
           data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
reg10 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37 + state_post37:year_c |state_post37|0|disease_year, 
            data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))

subset_state_data_1925_1943_long <- subset(subset_state_data_1925_1943_long, m_rate != 0)
reg11 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37|0|disease_year, 
            data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
reg12 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37 + state_post37:year_c |state_post37|0|disease_year, 
            data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))

table4_panel_b <- list(reg7, reg8, reg9, reg10, reg11, reg12)
latex_code <- texreg(table4_panel_b)
writeLines(latex_code, "table4_panel_b.tex")


##Panel C- State-level, excluding 1935 to 1937

state_data <- read_dta("data/datasets/state_data.dta")
# Filter the data to keep only the years between 1925 and 1943
state_data_1925_1943 <- state_data %>%
  filter((year >= 1925 & year <= 1934) | (year >= 1938 & year <= 1943))

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


# Run regression for each disease

reg13 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37|0|disease_year, 
             data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
reg14 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37 + state_post37:year_c |state_post37|0|disease_year, 
             data = filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
reg15 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37|0|disease_year, 
             data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
reg16 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37 + state_post37:year_c |state_post37|0|disease_year, 
              data = filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))

subset_state_data_1925_1943_long <- subset(subset_state_data_1925_1943_long, m_rate != 0)
reg17 <- felm(lnm_rate ~ treated:post37 + treated:year_c + treated + year_c + post37 |state_post37|0|disease_year, 
              data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
reg18 <- felm(lnm_rate ~ treated:post37 + treated:year_c:post37 + treated:year_c + treated + year_c + post37 + state_post37:year_c |state_post37|0|disease_year, 
              data = filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))


table4_panel_c <- list(reg13, reg14, reg15, reg16, reg17, reg18)
latex_code <- texreg(table4_panel_c)
writeLines(latex_code, "table4_panel_c.tex")




