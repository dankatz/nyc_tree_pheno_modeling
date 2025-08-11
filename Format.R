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
#batch_num <- 1 # Might need to extract multiple times to finish with a species

#### Pheno and XIS data ####
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

#Ends up having TreeID, species, Year, SOS_50, pixe_ID and pixels' lon and lat
tree_with_pixel <- tree_with_pixel %>%
  left_join(pixel_points, by = "pixel_ID") %>%  # Pixels' geometry here, not trees'
  arrange(pixel_ID) # In order to fit the metadata later!

#### XIS to Phenor ####

# Group by pixel_ID 
options(scipen = 999) # No scientific notation
pixel_groups <- tree_with_pixel  %>%
  group_by(pixel_ID) %>% # Still a data frame
  group_split() # Cut the data frame into several small data frames, and store them together as a list

# Split into smaller subset to extract and store, will be combined together later on
pixel_groups <- pixel_groups[1:2]

#### Extract Daymet Function (EPSG 4326) ####

extract_daymet_from_nc <- function(lon, lat, years, tile_dir, offset = 264, tile_ids = c(11753, 11754)) {
  vars <- c("vp", "prcp", "dayl")
  all_years <- sort(unique(c(years - 1, years))) # Extracting from previous year
  
  out_df <- purrr::map_dfr(
    all_years,  #Apply function(y) to each objects (i.e. year in all_years), will return a 365 rows tibble each time, being combined as a data frame
    function(y) {
      tibble(year = y,
             doy = 1:365,
             date = seq.Date(as.Date(paste0(y, "-01-01")), by = "day", length.out = 365))} #create a sequence from year-01-01, 365 days in total 
    )
  
  pt <- terra::vect(data.frame(x = lon, y = lat), geom = c("x", "y"), crs = "EPSG:4326") # Transferring  pixel_ID's lon&lat to SpatVector
  
  for (var in vars) {
    val_list <- list() # To store 365-day value for each year (E.g. val_list will be like 2017=c(XX, XX,...365 values), 2018=c(XX, XX,...365 values))
    
    for (yr in all_years) {
      val_year <- NA
      for (tile in tile_ids) {
        file <- file.path(tile_dir, paste0(var, "_", yr, "_", tile, ".nc")) #E.g. dayl_2018_11753.nc
        if (file.exists(file)) {
          r <- terra::rast(file) # Like a raster stack, each layer equals a day (365 layers)
          val <- terra::extract(r, pt)[1, -1] # From raster r, extract the value of "pt". Return a data frame, each columns represent a day. #Remove ID here (First col)
          val_year <- as.numeric(val)
          if (!all(is.na(val_year))) break
        }
      }
      val_list[[as.character(yr)]] <- val_year # Use yr as the name of list to store
    }
    
    var_mat <- do.call(cbind, val_list) #columns bind, combine each val_list (list) to a matrix 
    colnames(var_mat) <- paste0("year_", all_years)
    var_long <- as_tibble(var_mat) %>% # Convert a matrix to a tibble
      mutate(doy = 1:365) %>%
      pivot_longer(-doy, names_to = "year", values_to = var) %>% # wide format to long format (Column: doy, year, var)
      mutate(year = as.integer(gsub("year_", "", year)))
    
    # Merge to main df
    out_df <- left_join(out_df, var_long, by = c("year", "doy"))
  }
  return(out_df)
}

#### window function ####
get_365_window <- function(df, year, offset = 264) {
  start_date <- as.Date(sprintf("%s-01-01", year - 1)) + (offset - 1)
  length_needed <- 365

  if (lubridate::leap_year(year - 1)) { 
    length_needed <- 366}
  
  dates <- seq(start_date, length.out = length_needed, by = "day")
  
  df_window <- df %>% filter(date %in% dates) 
  
  return(df_window)
}

#### Format (Extract from Daymet, combine with XIS, and format to the Phenor-fitted list) ####
run_time <- system.time({

  site_data_list <- list()
  metadata_list <- list()  
  
  for(i in seq_along(pixel_groups)){
    cat(paste0("Pixel ", i, " ,"))
    group <- pixel_groups[[i]]
    px <- group$pixel_ID[1]
    lat <- group$lat[1]
    lon <- group$lon[1]
    years <- group$Year 
    sos   <- group$SOS_50
    
    # pixel's XIS
    climate_window <- xis_full %>%
      filter(pixel_ID == px,
             year %in% unique(c(years - 1, years))) 
    
    Ti_mat     <- purrr::map(years, ~ get_365_window(climate_window, .x)$temp_mean_C) %>% do.call(cbind, .)
    Tmini_mat  <- purrr::map(years, ~ get_365_window(climate_window, .x)$temp_min_C)  %>% do.call(cbind, .)
    Tmaxi_mat  <- purrr::map(years, ~ get_365_window(climate_window, .x)$temp_max_C)  %>% do.call(cbind, .)
    
    dm <- extract_daymet_from_nc(
      lon = lon,
      lat = lat,
      years = years,
      tile_dir = "C:/CapProject/NYC_daymet_reprj",  
      offset = 264)
    
    Pi_mat <- purrr::map(years, ~ get_365_window(dm, .x)$prcp) %>% do.call(cbind, .)
    vp_mat <- purrr::map(years, ~ get_365_window(dm, .x)$vp) %>% do.call(cbind, .)
    svp_mat <- 0.611 * exp((17.502 * Ti_mat) / (Ti_mat + 240.97))
    VPDi_mat <- svp_mat - (vp_mat / 1000)
    Li_mat <- purrr::map(years, ~ get_365_window(dm, .x)$dayl / 3600) %>% do.call(cbind, .)
    
    site_data_list[[i]] <- list(
      site = paste0("px_", px),
      location = c(lon, lat),
      doy = -102:262,
      ltm = rowMeans(Ti_mat, na.rm = TRUE),
      transition_dates = sos,
      year = years,
      Ti = Ti_mat,
      Tmini = Tmini_mat,
      Tmaxi = Tmaxi_mat,
      Li = Li_mat,
      Pi = Pi_mat,
      VPDi = VPDi_mat
    )
    # Store Tree_ID, Year, species
    metadata_list[[i]] <- group %>% select(Tree_ID, pixel_ID, Year, species, SOS_50)} #Pixel's lon lat
})





#### Combined all site running in this time/batch ####

combine_site_data <- function(site_data_list) {
  combined <- list(
    site     = "combined_site",
    location = site_data_list[[1]]$location,
    doy      = site_data_list[[1]]$doy,
    ltm      = site_data_list[[1]]$ltm
  )
  
  combined$transition_dates <- unlist(purrr::map(site_data_list, "transition_dates"))
  combined$year <- unlist(purrr::map(site_data_list, "year"))
  combined$Ti <- do.call(cbind, purrr::map(site_data_list, "Ti"))
  combined$Tmini <- do.call(cbind, purrr::map(site_data_list, "Tmini"))
  combined$Tmaxi <- do.call(cbind, purrr::map(site_data_list, "Tmaxi"))
  combined$Li <- do.call(cbind, purrr::map(site_data_list, "Li"))
  combined$Pi <- do.call(cbind, purrr::map(site_data_list, "Pi"))
  combined$VPDi <- do.call(cbind, purrr::map(site_data_list, "VPDi"))
  
  return(combined)
}

combined_site <- combine_site_data(site_data_list)

# Store combined site as RDS file
saveRDS(
  combined_site,
  file = file.path("C:/CapProject/25Summer/Extracted",
                   focal_species_under,
                   paste0(focal_species_under, batch_num, "_combined_site.rds"))
)


