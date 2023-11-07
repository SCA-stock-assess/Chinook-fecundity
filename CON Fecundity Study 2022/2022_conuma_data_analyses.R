# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "ggpmisc", "viridis", "RColorBrewer", "ggrepel", "readxl")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw(base_size = 20))
library(RColorBrewer)
library(viridis)
library(ggpmisc)
library(ggrepel)
library(readxl)


# Load and tidy data ------------------------------------------------------

# Data load and cleaning
fec <- read_xlsx("ConumaChinookFecundityData_2022.xlsx", 
                 skip = 2) |> # Skip the first two rows
  rename_with(~ str_replace(.x, "Sample", "Subsample") |> 
                str_replace("(?<=Subsample \\d) ", "_")) |> 
  pivot_longer(contains("_"),
               names_to = c("Subsample", "Measurement"),
               names_sep = "_",
               names_transform = str_trim) |> 
  select(-1, -`Estimated No Eggs`, -matches("Average|Diff"), -Comment) |> 
  filter(!is.na(Date),
         !grepl("Avg", Measurement)) |> 
  # Trim unit specifications from column names, abbreviate, and remove spaces
  rename_with(~ str_remove_all(.x, "\\(.*\\)") |> 
                abbreviate(minlength = 15) |> 
                str_replace_all("\\s", "_") |> 
                str_remove_all("_#")) |> 
  group_by(F_ID, Subsample) |> 
  mutate(Measurement = case_when(
    Measurement == "Weight" ~ "Egg Weight",
    grepl("^\\d\\d", Measurement) ~ gsub("\\d", "", Measurement),
    T ~ Measurement) |> 
      str_remove_all("\\(.*\\)") |> 
      abbreviate(minlength = 8),
    Subsample = str_replace(Subsample, "Subsample ", "ss"),
    ssEggCount = value[Measurement == "EggCount"],
    value = case_when(
      grepl("Lngth", Measurement) ~ value/10, 
      grepl("Weght", Measurement) ~ value/ssEggCount)) |>
  filter(!Measurement == "EggCount") |> 
  select(-ssEggCount) |> 
  pivot_wider(names_from = Subsample,
              values_from = value) |> 
  mutate(`ss%diff` = if_else(is.na(ss3),
                             abs(ss1 - ss2) / ss2,
                             (abs(ss1 - ss2) + abs(ss2 - ss3)) / ss2),
         est_ttl_egg_wt = TotalEggWeight / (`Est%ofEggRetntn`/100)) |> 
  rowwise() |> 
  mutate(ss_mean = mean(c_across(matches("ss\\d$")), na.rm = TRUE)) |> 
  group_by(F_ID) |> 
  mutate(est_fec = est_ttl_egg_wt / ss_mean[Measurement == "EggWeght"],
         mean_egg_di = ss_mean[Measurement == "EggLngth"]) |> 
  ungroup() %>%
  mutate(fec_resid = residuals(lm(est_fec ~ POH_Length)),
         di_resid = residuals(lm(mean_egg_di ~ POH_Length)),
         # Add columns containing lower and upper 95% CIs for the residuals 
         across(contains("resid"), 
                .fns = list(lwr = ~ quantile(.x, 0.025),
                            upr = ~ quantile(.x, 0.975)),
                .names = "{.fn}_{.col}")) 



# Plots -------------------------------------------------------------------

# Filtered datasets with outliers removed
fec_filt <- fec %>%
  filter(!(Measurement == "EggLngth" &
             (di_resid < lwr_di_resid |
             di_resid > upr_di_resid)),
         !(Measurement == "EggWeght" &
             (fec_resid < lwr_fec_resid |
                fec_resid > upr_fec_resid)))

# Egg diameter versus female length and fecundity
list(fec, fec_filt) |> 
  map_dfr(~ filter(.x, Measurement == "EggLngth"),
          .id = "id") |> 
  mutate(id = recode_factor(id, `1` = "Full", `2` = "Outliers\ntrimmed")) |> 
  ggplot(aes(POH_Length, ss_mean, colour = id)) +
  geom_point(data = filter(fec, Measurement == "EggLngth"),
             aes(fill = `ss%diff`, size = est_fec),
             colour = "black",
             shape = 21) +
  geom_smooth(method = "lm") +
  geom_text_repel(data = filter(fec, Measurement == "EggLngth"),
                  aes(label = ifelse(!between(di_resid, max(lwr_di_resid), max(upr_di_resid)),
                                     F_ID,
                                     "")),
                  colour = "black") +
  stat_poly_eq(use_label(c("eq", "R2", "n"))) +
  scale_fill_viridis_c(option = "rocket", 
                       end = 0.8,
                       labels = scales::percent) +
  scale_colour_manual(values = brewer.pal(n=5,"RdBu")[c(1,5)]) +
  labs(x = "Female POH length (mm)",
       y = "Mean egg diameter (mm)",
       colour = "Dataset for\nregression\nfit",
       size = "Estimated fecundity",
       fill = "Difference\nbtw. sub-\nsamples")



# Fecundity versus female length and egg diameter
list(fec, fec_filt) |> 
  map_dfr(~ filter(.x, Measurement == "EggWeght"),
          .id = "id") |> 
  mutate(id = recode_factor(id, `1` = "Full", `2` = "Outliers\ntrimmed")) |> 
  ggplot(aes(POH_Length, est_fec, colour = id)) +
  geom_point(data = filter(fec, Measurement == "EggWeght"),
             aes(fill = `ss%diff`, size = mean_egg_di),
             colour = "black",
             shape = 21) +
  geom_smooth(method = "lm") +
  stat_poly_eq(use_label(c("eq", "R2", "n")), vjust = 4, family = "Courier") +
  geom_text_repel(data = filter(fec, Measurement == "EggWeght"),
                  aes(label = ifelse(!between(fec_resid, max(lwr_fec_resid), max(upr_fec_resid)),
                                     F_ID,
                                     "")),
                  colour = "black") +
  scale_fill_viridis_c(option = "rocket", 
                       end = 0.9, 
                       trans = "log1p",
                       labels = scales::percent) +
  scale_colour_manual(values = brewer.pal(n=5,"RdBu")[c(1,5)]) +
  labs(x = "Female POH length (mm)",
       y = "Estimated fecundity",
       colour = "Dataset for\nregression\nfit",
       size = "Mean egg\ndiamater",
       fill = "Difference\nbtw. sub-\nsamples")

