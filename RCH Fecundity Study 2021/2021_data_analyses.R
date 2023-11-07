# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse","ggbeeswarm","ggpmisc","ggridges",
          "magrittr","reshape2","viridis","emmeans")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw(base_size = 20))
library(viridis)
library(ggridges)
library(ggpmisc)
library(ggbeeswarm)
library(ggExtra)
library(magrittr)
library(reshape2)
library(emmeans)

# Load, summarize, and plot data ------------------------------------------

# Data load and cleaning
fec <- read.csv("2021_fec_dat_full.csv", na.strings = c("", "N/A")) %>% 
  # resolved age: use CWT age to replace scale age where available
  mutate(resolved_age = if_else(is.na(resolved_age_cwt), resolved_age_scale, resolved_age_cwt),
         age = str_extract(resolved_age, "^[[:digit:]]") %>% as.factor(),
         date = as.Date(date, format = "%d-%b-%y"),
         est_egg_ret = est_egg_ret %>% 
           str_replace("\\?","") %>% 
           as.numeric(),
         fec_adj = round(fec_est/(est_egg_ret/100),0),
         poh.bin = case_when(
           poh %in% 506:547 ~ "(506,547]",
           poh %in% 548:588 ~ "(548,588]",
           poh %in% 589:629 ~ "(589,629]",
           poh %in% 630:670 ~ "(630,670]",
           poh %in% 671:711 ~ "(671,711]",
           poh %in% 712:752 ~ "(712,752]",
           poh %in% 753:793 ~ "(753,793]",
           poh %in% 794:834 ~ "(794,834]",
           poh %in% 835:875 ~ "(835,875]",
           poh %in% 876:916 ~ "(876,916]") %>% 
           as.factor(),
         across(fec_adj, .fns = list(log = log, sqrt = sqrt))) %>% 
  filter(!(age =="3" & poh > 700)) %>%  # Remove two outlier points where females were likely mis-aged
  mutate(fec_adj = if_else(poh > 700 & fec_est < 3000, fec_est/0.75, fec_adj)) # and fix two that were likely partially spent


# Summarize fecundity per bin
bin.fec <- fec %>% 
  group_by(poh.bin) %>% 
  summarize(across(c(fec_est,fec_adj), .fns = list(median = median, sd = sd)))

# Summarize fecundity per age
(age.fec <- fec %>% 
  group_by(age) %>% 
  summarize(across(c(fec_est,fec_adj), .fns = list(median = median, sd = sd)))
)

# Pivot table of ages by bin
fec %>% 
  count(poh.bin, age) %>% 
  pivot_wider(names_from = age,
              values_from = n,
              values_fill = 0) %>% 
  rowwise() %>% 
  mutate(n = sum(c_across(`3`:`5`)))

# Plot per bin and age
fec %>% 
  rename(`POH bin (mm)` = poh.bin, `Total age` = age) %>% 
  pivot_longer(cols = c(`Total age`,`POH bin (mm)`)) %>% 
  ggplot(aes(fec_est, reorder(value, desc(value)))) +
  facet_grid(name~., space = "free", scales = "free_y", switch = "y") +
  geom_density_ridges(quantile_lines = TRUE, 
                      quantiles = 2, alpha = 0.7,
                      jittered_points = TRUE) +
  labs(x = "Estimated fecundity", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(strip.placement = "outside")


# Plot of estimated, adjusted fecundity (and transformations) by POH length
fec %>% 
  melt(measure.vars = c("fec_adj","fec_adj_log","fec_adj_sqrt")) %>% 
  mutate(variable = as.factor(variable)) %>% 
  ggplot(aes(poh, value)) +
  facet_wrap(~variable, scales = "free_y", strip.position = "left", ncol=1,
             labeller = as_labeller(c(
               fec_adj = "Fecundity",
               fec_adj_log = "ln(Fecundity)",
               fec_adj_sqrt = "sqrt(Fecundity)"
             ))) +
  geom_point(aes(size = avg_egg_len, colour = age), alpha = 0.6) +
  geom_smooth(method = "lm") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  scale_colour_viridis_d(option = "mako", direction = -1, end = 0.8) +
  labs(y = NULL,
       x = "Female POH length (mm)",
       size = "Avg egg\ndiameter\n(mm)",
       colour = "Age") +
  theme(strip.placement = "outside",
        strip.background = element_blank())

# Same as above, but without the transformed fecundity data
ggMarginal(
  ggplot(fec, aes(poh, fec_adj)) +
    geom_point(aes(size = avg_egg_len, colour = age), alpha = 0.5) +
    geom_smooth(method = "lm", colour = "black") +
    stat_poly_eq(formula = y~x, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE) + 
    scale_radius(breaks = c(7:9), range = c(3.5,4.5)) +
    scale_colour_viridis_d(option = "mako", direction = -1, end = 0.8) +
    guides(colour = guide_legend(title.position = "top", 
                                 label.position = "bottom", 
                                 direction = "horizontal",
                                 override.aes = list(size = 6)),
           size = guide_legend(title.position = "top", 
                               label.position = "bottom", 
                               direction = "horizontal")) +
    labs(y = "Estimated fecundity",
         x = "Female POH length (mm)",
         size = "Avg egg\ndiameter\n(mm)",
         colour = "Age") +
    theme(legend.position = c(0.05,0.85),
          legend.justification = c(0,1),
          legend.box = "horizontal",
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black",fill=scales::alpha("white",0.8))),
  type = "density", 
  groupColour = TRUE, groupFill = TRUE,
  margins = "y")

  
# Plot egg length versus fecundity
ggplot(fec, aes(avg_egg_len, fec_adj)) +
  geom_point(aes(size = poh), alpha = 0.6) +
  geom_smooth(method = "lm")
# Messy, but smaller females clearly produce fewer and smaller eggs. 

# Plot egg weight versus fecundity
ggplot(fec, aes(avg_egg_wt, fec_adj)) +
  geom_point(aes(size = poh), alpha = 0.6) +
  geom_smooth(method = "lm")
# Messy, but smaller females clearly produce fewer and smaller eggs.



# Sample size per bin and infilled age comps based on historical data
# used to estimate batch fecundities per age
fec %>% 
  count(poh.bin) %>% 
  left_join(read.csv("age-comp_LU.csv"),
            by = c("poh.bin" = "bin_range")) %>% 
  select(-X) %>% 
  mutate(across(X31:X51, ~round(.x*n,0), .names = "{.col}_n")) %>% 
  left_join(bin.fec %>% select(1:2)) %>% 
  mutate(across(X31_n:X51_n, ~round(.x*fec_est_median,0), .names = "{.col}.eggs")) %>% 
  pivot_longer(cols = matches("^X\\d\\d_"),
               names_to = c("age","var"),
               names_pattern = "X(.*)_(.*)",
               values_to = "count") %>% 
  pivot_wider(names_from = var, values_from = count,
              id_cols = c(age,poh.bin)) %>% 
  group_by(age) %>% 
  summarize(across(n:n.eggs, sum)) %>% 
  mutate(fec_est_batch = n.eggs/n)





# Modelling ---------------------------------------------------------------

# Regression model of fecundity and egg size on female size
summary(poh.lm <- lm(fec_adj ~ age, fec))
plot(poh.lm) # Looks fine

# Use emmeans to get estimated fecundity values by age
emmeans(poh.lm, ~age, infer = c(TRUE,FALSE))




# Analysis of fecundities in prior years based on 2021 data ---------------


