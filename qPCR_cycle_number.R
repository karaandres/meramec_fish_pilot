### This code analyzes the qPCR data to determine PCR cycle number
### Last updated 5.5.2023 by Kara Andres (akara@wustl.edu)

### Clear work environment and load packages
rm(list = ls())
library(dplyr)
library(ggplot2)

# load datasets & merge by sample/well name
pilot_study_metadata <- read.csv("pilot_study_metadata.csv", header=TRUE)
amplification_dat <- read.csv("16S_pilot_diluted_2023-04-19_amplification.csv", header=TRUE)
dat <- merge(pilot_study_metadata, amplification_dat, by.x="qpcr_well", by.y="Well.Position")

# calculate delta-delta-Rn, plot results, calculate max value (optimal cycle #)
dat2 <- as.data.frame(dat %>%
                       group_by(qpcr_well) %>%
                       arrange(Cycle) %>%
                       mutate(Delta.Delta.Rn = Delta.Rn-lag(Delta.Rn, default=first(Delta.Rn))))
ggplot(dat2, aes(Cycle, Delta.Delta.Rn)) + 
  geom_point() +
  facet_wrap(~qpcr_well, ncol = 12)

dat3 <- as.data.frame(dat2 %>% group_by(qpcr_well) %>%
  filter(Delta.Delta.Rn == max(Delta.Delta.Rn)) %>% # filter the data.frame to keep row where x is maximum
  select(Cycle) %>%
  arrange(qpcr_well))

dat3 <- merge(pilot_study_metadata[,c(3,9)], dat3, by="qpcr_well")
hist(dat3$Cycle)
mean(dat3$Cycle)
# write.csv(dat3, "optimal_cycle.csv")  

