library(tidyverse)
library(lmerTest)

load('trial_by_trial.RData')
marker$task <- 'marker'; pointing$task <- 'pointing'
df <- bind_rows(pointing, marker)
df$trial <- as.numeric(df$trial)
df$correct <- as.numeric(df$correct)
rm(marker, pointing)

smry <- df %>% group_by(task, trial) %>% summarise(mean = (mean(correct, na.rm = T) * 100))

mycolors <- c('#4CAECD','#ED9E09')

theme_set(theme_classic(18))
p <- ggplot(smry, aes(x = factor(trial), y = mean, group = task)) +
  geom_line(aes(color = task), size = 1) + 
  geom_point(aes(color = task), size = 2) +
  coord_cartesian(ylim = c(0,100)) +
  geom_hline(yintercept = 50, linetype = 'dashed') +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_blank(), legend.position = c(0.5,0.85), legend.direction = 'horizontal') +
  labs(x = 'Trial', y = 'Percent Correct')
p

#####
## First trial analysis
####

t1 <- filter(df, trial == 1)
t1.point <- filter(t1, task == 'pointing')
t1.marker <- filter(t1, task == 'marker')

with(t1.point, binom.test(x = c(sum(correct ==1), sum(correct==0))))
with(t1.marker, binom.test(x = c(sum(correct ==1), sum(correct==0))))

#####
## learning across trials
####
marker <- filter(df, task == 'marker')
pointing <- filter(df, task == 'pointing')

marker.mod <- lmer(correct~trial + (1|name), data = marker)
summary(marker.mod)

pointing.mod <- lmer(correct~trial + (1|name), data = pointing)
summary(pointing.mod)