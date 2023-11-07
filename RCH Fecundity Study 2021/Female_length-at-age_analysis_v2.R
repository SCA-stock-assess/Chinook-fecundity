# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse","ggbeeswarm","magrittr", "ggridges")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw(base_size = 20))
library(ggbeeswarm)
library(magrittr)
library(ggridges)

#Function to pass over steps in a pipeline but print outputs
pass_thru <- function(data,fun) {print(fun(data)); data}

# Load and plot data ------------------------------------------------------

#Raw data read
f.len0 <- read.csv("length-at-age_female_CN_v2.csv") %>% 
  rename(
    FORK.LENGTH.NF = FORK.LENGTH..NOSE.FORK.,
    POH.LENGTH = POST.ORBITAL.HYPURAL..POH.
  )

# Investigate levels present in collection site variables and age data
f.len <- f.len0 %>% 
  pass_thru(. %>% count(PROJECT, SAMPLE_SOURCE, STATAREA)) %>% 
  pass_thru(. %>% count(RESOLVED.AGE, AGE_GR)) %>% 
  filter(str_detect(PROJECT, "Robertson")) %>% 
  mutate(age = case_when(str_detect(AGE_GR, "M") & RESOLVED.AGE != "" ~ RESOLVED.AGE,
                         RESOLVED.AGE == "" & str_detect(AGE_GR, "M") ~ NA_character_,
                         TRUE ~ AGE_GR)) %>% 
  pass_thru(. %>% count(RESOLVED.AGE, AGE_GR, age)) %>% 
  filter(!(age %in% c("No Age", "S1",""))) %>% 
  pass_thru(. %>% count(age)) %>% # Looks good. Will need to trim out the sub 2s
  pass_thru(. %>% count(is.na(POH.LENGTH), is.na(FORK.LENGTH.NF))) %>% #POH is most common length
  filter(POH.LENGTH > 500, !POH.LENGTH > 925) # Get rid of outlier observations

# How many years represented?
f.len %>% count(YEAR)

# Plot
ggplot(f.len %>% filter(!(str_detect(age, "2"))),
       aes(age, POH.LENGTH)) +
  geom_quasirandom(alpha = .1)

# Create dataframe with sample sizes
samples <- f.len %>% 
  filter(!(str_detect(age, "2|6"))) %>% 
  count(age)

# Density plot
ggplot(f.len %>% filter(!(str_detect(age, "2|6"))),
       aes(POH.LENGTH)) +
  geom_density(aes(fill = age), alpha = 0.7) +
  labs(y = "Density", x = "Post Orbital to Hypural Length (mm)",
         fill = "GR Age") +
  annotate("text", x = c(615,700,770), y = c(0.0105,0.01,0.0093), 
           label = c("N = 1,398","N = 9,382","N = 3,420"), colour = scales::hue_pal()(3)) +
  theme(legend.position = c(0.1,0.85),
        legend.background = element_rect(fill = "white",colour = "black"))
  

#Summary table
f.len %>% 
  group_by(age) %>% 
  summarize(across(c(POH.LENGTH, FORK.LENGTH.NF),median,na.rm=TRUE))

# Time series plot
ggplot(f.len, aes(YEAR, POH.LENGTH)) +
  geom_quasirandom(alpha = 0.2) +
  geom_smooth(method = "lm")


# Age distributions within length bins ------------------------------------

# Add bins
f.len %<>% mutate(POH.bin8 = cut(POH.LENGTH,8),
                  POH.bin10 = cut(POH.LENGTH,10),
                  POH.bin15 = cut(POH.LENGTH,15),
                  POH.bin20 = cut(POH.LENGTH,20),
                  POH.bin.used = case_when(
                    POH.LENGTH %in% 506:547 ~ "(506,547]",
                    POH.LENGTH %in% 548:588 ~ "(548,588]",
                    POH.LENGTH %in% 589:629 ~ "(589,629]",
                    POH.LENGTH %in% 630:670 ~ "(630,670]",
                    POH.LENGTH %in% 671:711 ~ "(671,711]",
                    POH.LENGTH %in% 712:752 ~ "(712,752]",
                    POH.LENGTH %in% 753:793 ~ "(753,793]",
                    POH.LENGTH %in% 794:834 ~ "(794,834]",
                    POH.LENGTH %in% 835:875 ~ "(835,875]",
                    POH.LENGTH %in% 876:916 ~ "(876,916]",
                    TRUE ~ NA_character_) %>% 
                    as.factor())


# Plot distribution within bins over time
YR.hist <- f.len %>% 
  filter(!str_detect(age, "2|6")) %>% 
  mutate(age = factor(age, levels = c("51","41","31"))) %>% 
  group_by(YEAR,POH.bin10,age) %>% 
  add_count() %>% 
  group_by(YEAR,age) %>% 
  mutate(N = sum(n),
         density = if_else(is.nan(n/N), 0, as.double(n/N))) %>% 
  ungroup()


ggplot(YR.hist, aes(y = POH.bin10, x = density)) +
  facet_grid(age~YEAR, switch = "x") +
  geom_bar(stat = "identity", aes(fill = n)) +
  geom_hline(yintercept = seq(0.5,10.5, by = 1), colour = "grey70") +
  theme_classic(base_size = 18) +
  labs(x = NULL, y = "POH length bin (mm)") +
  scale_fill_viridis_c(option = "inferno", end = 0.85, direction = -1) +
  #guides(fill = FALSE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0, "lines"),
        panel.background = element_rect(colour = "grey92"),
        strip.text = element_text(angle = 45, vjust = 1))


# Same plot, but as a heatmap
# First, make df with all counts by bin and year
YR.bins <- with(f.len %>% 
                  filter(!str_detect(age, "2|6")) %>% 
                  distinct(YEAR,POH.bin10, age), 
                expand.grid(
  YEAR = seq.int(min(YEAR),max(YEAR)),
  POH.bin10 = levels(POH.bin10),
  age = levels(as.factor(age))
)) %>% 
  left_join(f.len %>% 
              filter(!str_detect(age, "2|6")) %>% 
              group_by(YEAR, POH.bin10, age) %>% 
              count()) %>% 
  mutate(n = if_else(is.na(n), 0, as.double(n)),
         age = factor(age, levels = c("51","41","31"))) %>% 
  group_by(YEAR,age) %>% 
  mutate(N = sum(n),
         density = if_else(is.nan(n/N), 0, as.double(n/N)))

# Then make the heatmap plot
ggplot(YR.bins,aes(as.factor(YEAR),POH.bin10)) +
  geom_tile(aes(fill = density)) +
  geom_text(data = YR.bins %>% mutate(n = if_else(n == 0, "", as.character(n))), 
            aes(label = n), colour = "white") +
  #geom_hline(yintercept = seq(0.5,10.5, by = 1), 
  #           colour = "grey95", size = 0.2, lty = 3) +
  coord_cartesian(expand = FALSE) +
  scale_fill_viridis_c(option = "inferno", end = 0.9,
                       labels = scales::percent) +
  labs(x = NULL, y = "POH length bin (mm)") +
  facet_wrap(~age, ncol = 1, strip.position = "right",
             labeller = as_labeller(c(
               "51" = expression("5"[1]),
               "41" = expression("4"[1]),
               "31" = expression("3"[1])
             ))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_rect(colour = "black", fill = "white"))



# Calculate number of females of each age per bin
POH.comp <- f.len %>% 
  filter(!is.na(POH.LENGTH),!str_detect(age, "2|6")) %>% #Ignore/remove extreme ages
  pivot_longer(starts_with("POH.bin"),
               names_to = "bin_count",
               values_to = "bin_range",
               names_prefix = "POH.bin") %>% 
  group_by(bin_count, bin_range, age) %>% 
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n), bin_count = as.numeric(bin_count))

#write.csv(POH.comp, "bin_compositions.csv")

# Plot
POH.comp %>% 
  filter(is.na(bin_count), !is.na(bin_range)) %>% 
ggplot(aes(bin_range, freq, fill = age)) +
  #facet_wrap(~bin_count, scales = "free") +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = .5), size = 3) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "Binned post-orbital to hypural length (mm)", 
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save POH comps from 2021 chosen bins as a lookup table
(POH.comp 
  %>% filter(is.na(bin_count), !is.na(bin_range)) # select only the bins actually used in 2021
  %>% pivot_wider(names_from = age,
              values_from = freq, 
              id_cols = bin_range,
              values_fill = 0)
  %>% write.csv("age-comp_LU.csv")
)


# Ridgeline plot

# Calculate baseline period mean lengths
baseline <- f.len %>% 
  filter(!is.na(POH.LENGTH),!str_detect(age, "2|6"), YEAR %in% 1998:2003) %>%
  group_by(age) %>% 
  summarize(base_mean = mean(POH.LENGTH))
  
# Subset data and plot
f.len %>% 
  filter(!is.na(POH.LENGTH),!str_detect(age, "2|6")) %>%
  group_by(YEAR,age) %>% 
  add_count() %>% 
  group_by(age) %>% 
  ggplot(aes(x = POH.LENGTH, y = YEAR, fill = ..x.., group = YEAR)) +
  geom_density_ridges_gradient(bandwidth = 20, 
                               scale = 2,
                               rel_min_height = 0.01,
                               alpha = 0.4,
                               quantile_lines = TRUE,
                               quantile_fun = function(x,...)mean(x),
                               jittered_points = TRUE, 
                               point_alpha = 0.2,
                               point_size = 0.5) +
  geom_vline(data = baseline, aes(xintercept = base_mean),
             colour = "red", lty = 2) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  scale_y_reverse(breaks = seq.int(min(f.len$YEAR),max(f.len$YEAR)),
                  minor_breaks = NULL) +
  guides(fill = FALSE) +
  facet_wrap(~age, nrow = 1) +
  labs(x = "RCH Female Chinook POH length (mm)", y = NULL) +
  theme_minimal() +
  theme(panel.spacing = unit(0.5, "lines"),
        strip.background = element_rect(colour = "black", fill = "white"))

