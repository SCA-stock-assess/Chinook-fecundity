# Packages and functions --------------------------------------------------

pkgs <- c("tidyverse", "ggridges", "readxl", "writexl", "here", "httr", "jsonlite", "magrittr",
          "ggExtra", "ggpmisc")
#install.packages(pkgs)
#remotes::install_git("https://github.com/Pacific-salmon-assess/saaWeb")

library(tidyverse); theme_set(theme_bw(base_size = 18))
library(saaWeb)
library(ggridges)
library(ggpmisc)
library(ggExtra)
library(readxl)
library(writexl)
library(here)
library(httr)
library(jsonlite)
library(magrittr)

# Load age data from PADS online database ---------------------------------

# Age data from web services
# Thanks to Mike O'Brien for the code

# Declare function for extracting data from a single batch


# Extract web data 
sc_age_batches <- saaWeb::getAgeBatchList() |> 
  filter(
    SampleYear > 2020,
    Sector == "SC",
    Species == "Chinook",
    str_detect(ProjectName, "(?i)robertson|conuma")
  ) |> 
  pull(Id)


# Get ages
sc_ages <- saaWeb::getAgeBatchScaleResults(sc_age_batches) |> 
  janitor::clean_names()



# Load collection data from all study years -------------------------------


# Table of data files with sheet names
data_files <- list.files(
  here(),
  recursive = TRUE,
  full.names = TRUE
) |> 
  str_subset("~", negate = TRUE) |> 
  str_subset("\\.xlsx") |> 
  str_subset("(?i)data") |> 
  purrr::set_names() |> 
  map(
    ~ excel_sheets(.x) |> 
      as_tibble_col(column_name = "sheet")
  ) |> 
  list_rbind(names_to = "path")


# Load and clean all data
fec_data0 <- data_files |> 
  filter(str_detect(sheet, "(?i)^data$")) |> # Extract sheets named "data"
  mutate( # load data from excel files into a new column called "data"
    data = map2(
      path, 
      sheet, 
      ~read_excel(
        .x, 
        .y, 
        skip = 2
        )
      )
    ) |> 
  pull(data) |> # Extract the "data" column as a list
  map( # Keep only the necessary columns and wrangle the names and empty rows
    ~ .x |> 
      select(matches("(?i)date|length|weight|count|scale|whatman|id|retention|comment|clip")) |>
      select(!matches("(?i)avg|average")) |> # Remove averages (can calculate as/when needed)
      rename_with(
        ~ .x |> 
          tolower() |> 
          str_remove_all("[[:punct:]]") |> 
          str_trim() |> 
          str_replace_all("[[:space:]]+", "_") |> 
          str_replace_all("subsample", "sample")
      ) |> 
      filter(!is.na(date)) |> 
      mutate(year = format(date, "%Y"),
             f_id = paste0(year, "_", f_id),
             across(matches("length|weight|count|scale_book|whatman|retention"), as.numeric))
  ) |> 
  list_rbind() |> 
  pivot_longer( # Put the samples for egg weight and length into separate rows
    matches("sample_[[:digit:]]"),
    names_to = c("sample", "measure"),
    names_pattern = "sample_([[:digit:]])_(.*)"
  ) |> 
  pivot_wider(
    names_from = measure,
    values_from = value
  ) |> 
  mutate(
    site = case_when(
      year %in% c(2021, 2023) ~ "RCH",
      year %in% c(2022) ~ "CON",
      T ~ "Review"
    ),
    egg_length = `10_egg_length_mm` / 10,
    avg_f_egg_weight = sum(weight, na.rm = TRUE) / sum(egg_count, na.rm = TRUE),
    est_fecundity = as.integer(total_egg_weight_gm / avg_f_egg_weight),
    .by = f_id
  ) |> 
  rename("sample_weight" = "weight")


# Merge age data with original collection data
fec_data <- fec_data0 |> 
  mutate(across(everything(), as.character)) |> 
  left_join(
    select(sc_ages, matches("(?i)container|number|age")),
    by = c(
      "scale_book_no" = "container_id",
      "scale_book_fish_no" = "fish_number"
    )
  ) |> 
  mutate(across(everything(), parse_guess)) |> 
  rename_with(tolower) |> # Convert new columns to lowercase
  mutate(
    resolved_age_gr = case_when( 
      # Use CWT ages in preference to scale ages
      !is.na(resolved_scale_age_jt) ~ as.character(resolved_scale_age_jt),
      str_detect(gr_age, "[[:digit:]]{1}(?=M)") ~ paste0(as.numeric(str_extract(gr_age, "^[[:digit:]]")) + 1, 1),
      str_detect(gr_age, "[[:alpha:]]") ~ NA_character_,
      TRUE ~ gr_age
    ),
    comments = case_when(
      str_detect(comment, "^[[:digit:]]+$") ~ comments,
      is.na(comment) & is.na(comments) ~ NA_character_,
      TRUE ~ paste(comments, ";", comment)
    ),
    avg_f_egg_weight = as.numeric(avg_f_egg_weight)
  ) |> 
  select(-contains("cwt_ageresol"), -comment, -`10_egg_length_mm`)

# Add CWT data from EPRO?


# Collapse data to average over subsamples
fec_data_sum  <- fec_data |> 
  mutate(est_fec_corr = est_fecundity / (est_of_egg_retention/100)) |> # Correct for egg loss
  summarise(
    avg_egg_length = mean(egg_length, na.rm = TRUE),
    .by = c(
      date, 
      f_id,
      poh_length_mm,
      ad_clip_yn,
      comments,
      year,
      resolved_age_gr,
      site,
      avg_f_egg_weight,
      est_fecundity,
      est_fec_corr
    )
  )

# Save full and summarized data into one excel file
list(fec_data, fec_data_sum) |> 
  map( # Fix up column names so they look nicer in Excel
    ~ .x |> 
      rename_with(
        ~ .x |> 
          str_replace_all("_", " ") |> 
          str_to_title() |> 
          str_replace_all("(?i)poh|yn|gr|eu|id", str_to_upper) |> # Capitalize all the acronyms
          str_replace_all("(?i)mm", str_to_lower) # Convert millimeters (mm) to lowercase
      )
  ) |> 
  set_names(c("Full data", "Summarized by F ID")) |> # Specify sheet names
  write_xlsx(
    path = here(
      paste0(
        "R-OUT_Fecundity-Data_all-years_", 
        Sys.Date(), 
        ".xlsx"
      )
    )
  )


# Plot fecundity data -----------------------------------------------------


# Holistic plot of stocks and years
(p1 <- fec_data_sum |> 
   mutate(avg_f_egg_weight = avg_f_egg_weight * 1000) |> # Convert to mg
   pivot_longer(c(est_fec_corr, avg_f_egg_weight, avg_egg_length)) |> 
   ggplot(aes(
     x = poh_length_mm, 
     y = value, 
     colour = site)
   ) +
   facet_wrap(
     ~name, 
     ncol = 1, 
     strip.position = "left",
     scales = "free_y",
     labeller = as_labeller(
       c(
         "avg_egg_length" = "Mean egg diameter (mm)",
         "avg_f_egg_weight" = "Mean egg weight (mg)",
         "est_fec_corr" = "Estimated fecundity"
       )
     )
   ) +
   geom_point(
     aes(shape = as.factor(year)),
     alpha = 0.5
   ) +
   geom_smooth(method = "lm") +
   stat_poly_eq(
     formula = y ~ x,
     aes(label = ..rr.label..),
     parse = TRUE,
     label.y = c(0.85, 0.95)) +
   guides(shape = guide_legend(override.aes= list(size = 3))) +
   labs(
     y = NULL,
     x = "Post-orbital to hypural length (mm)",
     shape = "year"
   ) +
   theme(
     strip.placement = "outside",
     strip.background.y = element_blank()
   )
)


# Ridgeline plot showing fecundity by age
(p2 <- fec_data_sum |> 
    filter(!is.na(resolved_age_gr)) |> 
    ggplot(aes(est_fec_corr, site)) +
    facet_wrap(
      ~resolved_age_gr,
      ncol = 1,
      strip.position = "right",
      labeller = as_labeller(function(x) (paste("Age", str_extract(x, "[[:digit:]]"))))
    ) +
    geom_density_ridges(
      aes(point_colour = as.character(year)),
      jittered_points = TRUE,
      position = "raincloud",
      alpha = 0,
      colour = NA,
      point_alpha = 0.7
    ) +
    geom_density_ridges(
      alpha = 0.7,
      quantile_lines = TRUE,
      quantiles = 2,
      vline_colour = "red"
    ) +
    scale_discrete_manual(
      aesthetics = "point_colour",
      values = RColorBrewer::brewer.pal(length(unique(fec_data_sum$year)), "Dark2"),
      name = "year",
      guide = guide_legend(override.aes = list(point_size = 3))
    ) +
    labs(
      x = "Estimated fecundity",
      y = "Hatchery stock sampled"
    ) +
    theme(
      legend.position = c(0.98,0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(colour = "black")
    )
)


# Plots showing fecundity by length and age
(p3 <- fec_data_sum |> 
    filter(!is.na(resolved_age_gr),
           !(poh_length_mm > 700 & resolved_age_gr == "31")) |> # Remove some outliers
    mutate(
      group = paste0(site, year),
      resolved_age_gr = str_extract(resolved_age_gr, "^[[:digit:]]")
    ) %>%
    split(.$group) |> 
    map(
      function(data) {
        title <- paste(unique(data$site), unique(data$year))
        
        plot <- data |> 
          ggplot(aes(poh_length_mm, est_fec_corr)) +
          geom_point(aes(colour = resolved_age_gr, size = avg_egg_length),
                     alpha = 0.75) +
          geom_smooth(method = "lm", colour = "black") +
          stat_poly_eq(formula = y~x, 
                       aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                       parse = TRUE) + 
          scale_radius(
            breaks = c(6:9), 
            range = c(3, 4.5),
            guide = guide_legend(
              title.position = "top",
              label.position = "bottom",
              label.vjust = 3
            )
          ) +
          scale_colour_viridis_d(
            option = "mako", 
            direction = -1, 
            end = 0.8,
            guide = guide_legend(
              title.position = "top",
              label.position = "bottom",
              override.aes = list(size = 4),
              label.vjust = 3
            )
          ) +
          labs(
            x = "Post-orbital to hypural length (mm)",
            y = "Estimated fecundity",
            size = "Mean egg diameter (mm)",
            colour = "Total age (years)",
            #title = title
          ) +
          theme(legend.position = "bottom")
        
        marginal_plot <- ggMarginal(
          plot,
          type = "density", 
          groupColour = TRUE, 
          groupFill = TRUE,
          margins = "y"
        )
        
        return(marginal_plot)
      }
    )
)


# Save plots
# Demographic correlations
ggsave(
  plot = p1,
  filename = here(
    "plots",
    "R-PLOT_fecundity_egg-size_length_correlations.png"
  ),
  height = 8,
  width = 7,
  units = "in"
)


# Ridgelines
ggsave(
  plot = p2,
  filename = here(
    "plots",
    "R-PLOT_fecundity_by_age_ridges.png"
  ),
  height = 8,
  width = 8,
  units = "in"
)


# Marginal histogram plots
p3 |> 
  iwalk(
    ~ ggsave(
      plot = .x,
      filename = paste0(
        here("plots"),
        "/R-PLOT_fecundity_by_length_",
        .y,
        ".png"
      ),
      width = 8,
      height = 6,
      units = "in"
    )
  )
