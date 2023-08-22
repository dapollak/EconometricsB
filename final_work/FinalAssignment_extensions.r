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


## mortality trends

national_data <- read_dta("national_data.dta")
subset_national_data <- national_data[, c("mmr", "influenza_pneumonia_total", "scarlet_fever_tot", "tuberculosis_total", "year")]

#Reshape long
subset_national_data_long <- subset_national_data %>%
  pivot_longer(cols = c(mmr, influenza_pneumonia_total, scarlet_fever_tot, tuberculosis_total),
               names_to = "disease",
               values_to = "m_rate") 

# Create new variables
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


##parallel trend estimation:

state_data <- read_dta("state_data.dta")
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
  mutate(treated = as.numeric(disease == "mmr" | disease == "infl_pneumonia_rate" | disease == "scarfever_rate"),
         year_c = year - 1937,
         lnm_rate = log(m_rate))

subset_state_data_1925_1943_long$disease_year <- as.numeric(interaction(subset_state_data_1925_1943_long$disease, subset_state_data_1925_1943_long$year))  

# Run regression for each disease- omitted 1925
ev_study_mmr <- feols(lnm_rate~ i(year,treated,ref=1925) | year + treated + state^year, cluster="disease_year",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
plt1 = coefplot(ev_study_mmr, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.2,1.2))
summary(ev_study_mmr)

ev_study_inf_pne <- feols(lnm_rate~ i(year,treated,ref=1925) | year + treated + state^year, cluster="disease_year",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
plt2 = coefplot(ev_study_inf_pne, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.2,1.2))
summary(ev_study_inf_pne)

ev_study_scarlet <- feols(lnm_rate~ i(year,treated,ref=1925) | year + treated + state^year, cluster="disease_year",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
plt3 = coefplot(ev_study_scarlet, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.2,1.2))
summary(ev_study_scarlet)


# Run regression for each disease- omitted 1937
reg19 <- feols(lnm_rate~ i(year,treated,ref=1937) | year + treated + state^year, cluster="disease_year",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
plt1 = coefplot(reg19, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.2,1.2),title = "Maternal Mortality")

reg20 <- feols(lnm_rate~ i(year,treated,ref=1937) | year + treated + state^year, cluster="disease_year",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
plt2 = coefplot(reg20, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.2,1.2),title = "Influenza and Pneumonia")

reg21 <- feols(lnm_rate~ i(year,treated,ref=1937) | year + treated + state^year, cluster="disease_year",
               data=filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
plt3 = coefplot(reg21, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-1.2,1.2),title = "Scarlet Fever")



# Run regression for each disease but with level instead of log- omitted 1925
ev_study_mmr_level <- feols(m_rate~ i(year,treated,ref=1925) | year + treated + state^year, cluster="disease_year",
                      data=filter(subset_state_data_1925_1943_long, disease %in% c("mmr", "tb_rate")))
summary(ev_study_mmr_level)
plt1 = coefplot(ev_study_mmr_level, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-500,10))

ev_study_inf_pne_level <- feols(m_rate~ i(year,treated,ref=1925) | year + treated + state^year, cluster="disease_year",
                          data=filter(subset_state_data_1925_1943_long, disease %in% c("infl_pneumonia_rate", "tb_rate")))
summary(ev_study_inf_pne_level)
plt2 = coefplot(ev_study_inf_pne_level, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-0.1,0.1))

ev_study_scarlet_level <- feols(m_rate~ i(year,treated,ref=1925) | year + treated + state^year, cluster="disease_year",
                          data=filter(subset_state_data_1925_1943_long, disease %in% c("scarfever_rate", "tb_rate")))
summary(ev_study_scarlet_level)
plt3 = coefplot(ev_study_scarlet_level, outerCI=1.96, innerCI=1.645, horizontal = TRUE, ylim=c(-0.1,0.1))





