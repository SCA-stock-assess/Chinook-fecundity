# Packages & functions ----------------------------------------------------

pkgs <- c("VCA","VFP","simr","tidyverse","splitstackshape")
install.packages(pkgs)

library(VCA)
library(VFP)
library(simr)
library(tidyverse); theme_set(theme_bw(base_size=20))
library(splitstackshape)

set.seed(3)

# Create dataset with biostandard assumptions by age ----------------------


# CV for RCH Chinook estimated below based on data from Healey & Heard (1984)
(RCH.CV = ((6799-2943)/4)/4452) # Authors provide only the mean (4452) & range (2943-6799)
# Uses range rule to calculate SD (SD = range/4), which only works for gaussian distributions.

# Create a dataframe with the group-level parameters
params <- data.frame(
  age = c(3,4,5), 
  fecundity = c(2500,3700,4200),
  # CV for RCH Chinook estimated below based on Kaufman et al. (2009) data on Mokelumne river CN
  cv_fec = c(1348/4185, 1205/5838, 784/5994) # appears to decrease with age
) %>% 
  mutate(sd_fec = cv_fec*fecundity, 
         age = as.factor(age))

# Use the parameters to infill a dataframe of simulated fecundity values
sim.dat <- map_df(1:50, ~params) %>% 
  mutate(est_fec = map2(fecundity, sd_fec, ~rnorm(1,.x,.y)),
         across(everything(), unlist))

# Add mean and sd from sim.dat to the source parameters
comp <- sim.dat %>% 
  group_by(age) %>% 
  summarize(fecundity = mean(est_fec), sd_fec = sd(est_fec)) %>% 
  mutate(source = "simulation") %>% 
  bind_rows(params %>% 
              select(age, fecundity, sd_fec) %>% 
              mutate(source = "parameters"))


# Plot of expected and simulated data distribution
(pred.p <- ggplot(comp, aes(age, fecundity, colour = source)) +
    geom_pointrange(aes(ymin = fecundity - (sd_fec)*1.96,
                        ymax = fecundity + (sd_fec)*1.96),
                    size = 1.5, position = position_dodge(width = .3))
)
# Perfect

# simr power analysis --------------------------------------------------------------

# First, fit a model to the simulated data
mod <- lm(est_fec ~ age, sim.dat)
summary(mod)

# Add levels
mod2 <- extend(mod, along = "est_fec", n = 30)

# Show power curve by sample size
plot(powerCurve(mod2, along = "est_fec", breaks = seq(2,30, by = 2)), 
     xlab = "Sample size per age class")



# Precision curves based on manual simulation -----------------------------

# Create dataframe with groups at varying sizes N per age
n.range <- map_dfr(2:70,~expandRows(params, count = .x, count.is.col= FALSE), .id = "N") %>% 
  mutate(est_fec = map2(fecundity, sd_fec, ~rnorm(1,.x,.y)),
         across(everything(), unlist),
         N = (as.integer(N)+1) %>% as.factor)

# Make sure it worked...
n.range %>% count(N) #Success!

# Baby example with single model run on subset of the data
bb.mod <- lm(est_fec~age, data = subset(n.range, N == "15")) # Run the model
confint(bb.mod) %>% # Mess around with the 'confint' output to get it in a nice format
  as_tibble() %>% 
  mutate(width = `97.5 %` - `2.5 %`) %>% 
  select(width) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(age3 = V1, age4 = V2, age5 = V3)


# Summarize and plot results from all the models
n.range %>% 
  split(.$N) %>% 
  map(~lm(est_fec~age, data=.x)) %>% 
  map_dfr(~confint(.) %>%
            as_tibble() %>% 
            mutate(width = `97.5 %` - `2.5 %`) %>% 
            select(width) %>% 
            t() %>% 
            as.data.frame() %>% 
            rename(age3 = V1, age4 = V2, age5 = V3)
          ) %>%
  mutate(N = seq.int(2,70,by = 1)) %>% 
  pivot_longer(cols = age3:age5,
               names_to = "age", values_to = "CI_width") %>% 
  ggplot(aes(N, CI_width, colour = age)) +
  #geom_line(position = position_dodge(width = 1)) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"),
              position = position_dodge(width = 1), alpha = 0)
# These curves are based on just a single run of the simulation! Interpret with caution. 
