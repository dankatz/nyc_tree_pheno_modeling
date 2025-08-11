library(tidyverse)
library(daymetr)
library(dplyr)
library(lubridate)
library(arrow)
library(phenor)
library(sf)
library(purrr)
library(terra)
library(rbenchmark)
library(ggspatial)

setwd("C:/CapProject/25Summer")
focal_species <- "Quercus rubra"
focal_species_under <- "Quercus_rubra"


#### 1. Read in data, prepare "tree_with_pixel" ####
df <- read_csv("C:\\CapProject\\TreePheno\\tree_pheno_pointextract_polyid_all_output_2017_2024_copy.csv")

# Filter valid observations
# Tree observations
df <- df %>% 
  filter(species == focal_species, R2 > 0.8, tpstructur == "Full", tpconditio %in% c("Good", "Excellent")) %>%
  filter(dbh > 3.94) %>% #only take big trees i.e., those with dbh >10 cm
  filter(!Year %in% c(2017, 2024)) %>% 
  select(Point_ID, species, Lon, Lat, Year, SOS_50) %>%
  rename(Tree_ID = Point_ID) 

df_summary <- df %>% 
  group_by(species, Year) %>% 
  summarize(p5 = quantile(SOS_50, 0.05),
            p95 = quantile(SOS_50, 0.95))

df <- left_join(df, df_summary) %>% 
  filter(SOS_50 > p5 & SOS_50 < p95)

tree_sf <- st_as_sf(df, coords = c("Lon", "Lat"), crs = 4326)

xis_full <- arrow::read_parquet("xis_full_cleaned.parquet")


# Assign trees to nearest pixel

pixel_points <- xis_full %>% select(pixel_ID, lon, lat) %>% distinct()

xis_points_sf <- st_as_sf(pixel_points, coords = c("lon", "lat"), crs = 4326)

tree_with_pixel <- st_join(tree_sf, xis_points_sf,
                           join = st_nearest_feature,
                           left = TRUE) %>% st_drop_geometry() %>% select(Tree_ID, species, Year, SOS_50, pixel_ID)

#End up having Tree_ID, species, Year, SOS_50, pixe_ID and pixels' lon and lat
tree_with_pixel <- tree_with_pixel %>%
  left_join(pixel_points, by = "pixel_ID") %>%  # Pixels' geometry here, not trees'
  arrange(pixel_ID) # In order to fit the order of metadata later!

options(scipen = 999) # No scientific notation
pixel_groups <- tree_with_pixel  %>%
  group_by(pixel_ID) %>% # Still a data frame
  group_split()   # Cut the data frame into several small data frames, and store them together as a list
  
#### 2. Read in RDS filed created in "Format.R" ####
#Read in rds
rds_files <- list.files(path = paste0("C:\\CapProject\\25Summer\\Extracted\\", focal_species_under), 
                        pattern = paste0(focal_species_under, ".*_combined_site\\.rds$"), 
                        full.names = TRUE)

site_data_list <- map(rds_files, readRDS)


#### 3. Combine all batches and run models (whole NYC) ####
runtime_phenor <- system.time({
  combine_all_sites <- function(site_data_list) {
    combined <- list(
      site = focal_species_under,
      location = site_data_list[[1]]$location,
      doy = site_data_list[[1]]$doy,
      ltm = site_data_list[[1]]$ltm
    )
    
    combined$transition_dates <- unlist(map(site_data_list, "transition_dates"))
    combined$year <- unlist(map(site_data_list, "year"))
    combined$Ti <- do.call(cbind, map(site_data_list, "Ti"))
    combined$Tmini <- do.call(cbind, map(site_data_list, "Tmini"))
    combined$Tmaxi <- do.call(cbind, map(site_data_list, "Tmaxi"))
    combined$Li <- do.call(cbind, map(site_data_list, "Li"))
    combined$Pi <- do.call(cbind, map(site_data_list, "Pi"))
    combined$VPDi <- do.call(cbind, map(site_data_list, "VPDi"))
    
    return(combined)
  }
  
  all_sites <- combine_all_sites(site_data_list)
  
  models <- c("null", "LIN", "TT", "PTT", "M1", "AT", "SQ", "PA")
 
  results_all_sites <- list()
  
  for (model in models) {
    cat("Running", model, "\n")
    
    if (model == "null") {
      predicted <- null(all_sites)
      results_all_sites[[model]] <- list(
        predicted = predicted,
        measured = all_sites$transition_dates,
        rmse = sqrt(mean((predicted - all_sites$transition_dates)^2)),
        aic = list(AIC = NA, AICc = NA),
        model = "null"
      )
    } else {
      # Other models
      results_all_sites[[model]] <- pr_fit(
        model = model,
        data = all_sites,
        method = "GenSA",
        plot = FALSE
      )
    }
  }
})

# Store "results_all_sites" into "Model_Result" file, it's species result across whole NYC
saveRDS(results_all_sites,
        file = file.path("C:/CapProject/25Summer/Model_Result", paste0(focal_species_under, "_all_sites_results.rds")))

NYC_results <- readRDS(paste0("C:/CapProject/25Summer/Model_Result/", focal_species_under,"_all_sites_results.rds"))
                       
#### 4. Export result_summary and metadata (csv) ####
metadata_trees_NYC <- tree_with_pixel  %>% #metadata_trees_NYC is actually tree_with_pixel
  select(Tree_ID, pixel_ID, species, Year, lon, lat, SOS_50 ) %>%
  rename(transition_dates = SOS_50)

extract_summary <- function(results, site_name = focal_species_under) {
  tibble(
    Site = focal_species_under,
    Model = names(results),
    RMSE = purrr::map_dbl(results, "rmse"),
    AIC = purrr::map_dbl(results, ~ .$aic$AIC),
    AICc = purrr::map_dbl(results, ~ .$aic$AICc)
  )
}

results_summary <- extract_summary(results_all_sites) 

write.csv(results_summary,
  file = file.path("C:/CapProject/25Summer/Result_summary",
                   paste0(focal_species_under, "_model_summary.csv")),  row.names = FALSE)

df_all <- tibble(
  Tree_ID  = metadata_trees_NYC$Tree_ID,
  pixel_ID = metadata_trees_NYC$pixel_ID,
  Year     = metadata_trees_NYC$Year,
  Species  = metadata_trees_NYC$species,
  lon      = metadata_trees_NYC$lon,
  lat      = metadata_trees_NYC$lat,
  Observed = metadata_trees_NYC$transition_dates
)

for (model in names(results_all_sites)) {
  df_all[[model]] <- results_all_sites[[model]]$predicted
}

write.csv(df_all,
  file = file.path("C:/CapProject/25Summer/Model_Result",
                   paste0(focal_species_under, "_models_result.csv")),  row.names = FALSE)

#Make sure data frame merge  correctly
results_all_sites[["TT"]]$measured
metadata_trees_NYC$transition_dates
tree_with_pixel$SOS_50

