library(tidyverse)

data <- read.csv("/Users/sherront/Downloads/Safari Downloads/Sample Metadata.csv")
library(ggplot2)
theme_set(theme_bw())
ggplot(data, aes(x = case_submitter_id, y=pathologic_stage)) + 
  geom_point(stat='identity', aes(col=gender), size=3)  +
  theme(legend.position="none") + 
  scale_color_manual(name="gender", 
                     labels = c("male", "female"), 
                     values = c("male"="#0064d4", "female"="#ec4cf6")) +
  scale_x_discrete(limits=data$case_submitter_id)+
  ylab("Pathological Stage")+
  xlab("Patient Case ID")+
  coord_flip()