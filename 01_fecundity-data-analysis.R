# Packages and functions --------------------------------------------------

pkgs <- c("tidyverse", "readxl", "here", "httr", "jsonlite", "magrittr")
#install.packages(pkgs)

library(tidyverse); theme_set(theme_bw(base_size = 18))
library(readxl)
library(here)
library(httr)
library(jsonlite)
library(magrittr)

# Load age data from PADS online database ---------------------------------

# Age data from web services
# Thanks to Mike O'Brien for the code

# Declare function for extracting data from a single batch
get_age_detail <- function(batch_id){
  url <- paste0('http://pac-salmon.dfo-mpo.gc.ca/Api.CwtDataEntry.v2/api/AgeBatchDetail/ScaleAgeDetail/', 
                batch_id)
  x <- httr::GET(url, authenticate(':', ':','ntlm'))
  y <- batch_df_sc_sock <- jsonlite::fromJSON(httr::content(x, "text"))
  return(y)
}


# Extract web data 
sc_ages <- httr::GET(
  'http://pac-salmon.dfo-mpo.gc.ca/Api.CwtDataEntry.v2/api/Report/GetAgeBatchList', 
  httr::authenticate(':', ':','ntlm')
) %>%
  # Get list of available age result batches
  httr::content(x = .,  'text') |> 
  jsonlite::fromJSON() |>
  # Filter batch list to keep only Chinook samples from 2021 onward
  filter(
    Species  == 'Chinook',
    Sector == "SC",
    SampleYear > 2020,
    str_detect(ProjectName, "(?i)robertson|conuma")
  ) |> # Keep only Robertson and Conuma samples
  pull(Id) %>% # Extract column of age batch Ids from resulting dataframe
  purrr::map_dfr(., get_age_detail) |> # Run each batch ID through the function to extract its aging data
  mutate(across(matches("(?i)container|fishnumber"), as.numeric)) |> # Convert book and cell # to numeric
  rename(GR_Age = GrAge) 



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
  left_join(
    select(sc_ages, matches("(?i)container|number|age")),
    by = c(
      "scale_book_no" = "ContainerId",
      "scale_book_fish_no" = "FishNumber"
      )
  ) |> 
  rename_with(tolower) |> # Convert new columns to lowercase
  mutate(
    resolved_age_gr = case_when( # Use CWT ages in preference to scale ages
      !is.na(resolved_scale_age_jt) ~ as.character(resolved_scale_age_jt),
      TRUE ~ gr_age
    ) %>%
      if_else( # Convert partial ages to sub-1s
        str_detect(., "[[:digit:]]{1}(?=M)"),
        paste0(
          as.numeric(str_extract(., "^[[:digit:]]")) + 1, # Add 1 to marine age
          1
        ),
        .
      ),
    comments = case_when(
      str_detect(comment, "^[[:digit:]]+$") ~ comments,
      is.na(comment) & is.na(comments) ~ NA_character_,
      TRUE ~ paste(comments, ";", comment)
    )
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

# Plot fecundity data -----------------------------------------------------

# Holistic plot of stocks and years
fec_data_sum |> 
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
    aes(shape = year),
    alpha = 0.5
    ) +
  geom_smooth(method = "lm") +
  labs(y = NULL) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank())
