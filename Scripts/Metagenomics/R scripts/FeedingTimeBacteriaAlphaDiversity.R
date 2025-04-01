library(tidyverse)

diversitydata <- read_csv("FeedingTimeBacteriaAlphaDiversity.csv")



library(lme4)
library(lmerTest)  # Note, we might have to install this first! 


mem_rt_4 <- lmer(Diversity ~ TimePoint + Age + BMI + Antibiotic + Gender + (1|PatientID) + (1|SmokingStatus), data = diversitydata)
coef(summary(mem_rt_4))
head(coef(mem_rt_4)$PatientID)
summary(mem_rt_4)

t.test(Diversity ~ TimePoint, data=diversitydata, paired=TRUE)$statistic

dotplot(PatientID ~ resid(mem_rt_4), data = diversitydata, abline = list(v = 0))
summary(resid(mem_rt_4))
        
        
