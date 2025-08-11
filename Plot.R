library(tidyverse)
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

#### Read in the final CSV file ####
df_all <- read.csv(file = file.path("C:/CapProject/25Summer/Model_Result",
                   paste0(focal_species_under, "_models_result.csv"))
)


best_model <- "TT" #Remember switching to the best performed model

#### Scatter Plot (Predicted vs Observed)####
df_plot <- data.frame(
  Observed = df_all$Observed,
  Predicted = df_all[[best_model]] 
)

p <- ggplot(df_plot, aes(x = Predicted, y = Observed)) +
  #geom_point(alpha = 0.2, color = "black", size = 2) + 
  geom_hex(bins = 15) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "sqrt") +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.8)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1.3) +
  theme_light(base_size = 15) +
  labs(
    title = paste0("Predicted vs Observed - ", focal_species, " (", best_model, ")"),
    x = "Predicted DOY",
    y = "Observed DOY",
    fill = "Count"
  )+
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(30, 40, 30, 40)
  )

p

ggsave(filename = paste0(focal_species_under, "_Pred_Obs.png"),
  plot = p,  width = 12,  height = 10,  dpi = 600)

#### Scatter Plot Facet ####
df_alba <- read_csv("C:/CapProject/25Summer/Model_Result/Quercus_alba_models_result.csv") %>%
  mutate(Species = "Quercus alba",
         Predicted = TT)

df_robur <- read_csv("C:/CapProject/25Summer/Model_Result/Quercus_robur_models_result.csv") %>%
  mutate(Species = "Quercus robur",
         Predicted = AT)

df_rubra <- read_csv("C:/CapProject/25Summer/Model_Result/Quercus_rubra_models_result.csv") %>%
  mutate(Species = "Quercus rubra",
         Predicted = TT)

df_velutina <- read_csv("C:/CapProject/25Summer/Model_Result/Quercus_velutina_models_result.csv") %>%
  mutate(Species = "Quercus velutina",
         Predicted = AT)

# Combine those species to the same dataframe
df_all_species <- bind_rows(df_alba, df_robur, df_rubra, df_velutina)

df_plot <- df_all_species %>%
  select(Species, Observed, Predicted)

p_facet <- ggplot(df_plot, aes(x = Predicted, y = Observed)) +
  geom_hex(bins = 15) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "sqrt") +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1.3) +
  facet_wrap(~ Species, ncol = 2) + 
  labs(
    x = "Predicted (DOY)",
    y = "Observed (DOY)"
  )+
  theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid = element_line(color = "gray90"),
    panel.spacing = unit(1, "lines")
  )


ggsave(filename = "All_Species_Pred_Obs.png",
       plot = p_facet,  width = 12,  height = 10,  dpi = 900)

#### Prediction Map ####
df_all_sf <- df_all %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

#load in nyc boundary polygon
nyc_boundary_raw <- st_read("C:/CapProject/nyc_boundary_polygon/nybb.shp")

nyc_boundary <- nyc_boundary_raw %>%
  st_transform(crs = 4326) %>%
  st_union()  #combine the different boroughs

# NYC bounding box -> difference -> sf
nyc_boundary_box <- st_as_sf(st_as_sfc(st_bbox(nyc_boundary)),  crs = st_crs(nyc_boundary))
nyc_boundary_invert <- st_difference(st_geometry(nyc_boundary_box),st_geometry(nyc_boundary))

plot_map_facet <- function(df_sf, model_col = "AT") {
  data_model <- df_sf %>%
    filter(Species == focal_species)
  
  ggplot() +
    ggthemes::theme_few() +
    theme(
      panel.background = element_rect(fill = "gray94", color = NA),
      plot.title = element_text(hjust = 0.5)
    ) +
    
    geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +
    geom_sf(data = nyc_boundary, fill = "white", color = "black") +
    
    geom_sf(data = data_model, aes_string(color = model_col), size = 1.3, alpha = 0.8) +
    scale_color_gradientn(colors = c("red3", "yellow", "darkgreen"))+
    labs(
      title = paste(focal_species, "- Predicted (",model_col,")"),
       color = "DOY"
    )+
    annotation_scale(location = "br") +
    facet_wrap(~ Year)
}

pred_map <- plot_map_facet(df_all_sf, model_col = best_model)

pred_map

ggsave(
  filename = paste0(focal_species_under, "_pred_map.png"),
  plot = pred_map,  width = 12,  height = 8,  dpi = 600)


#### Residuals Map ####
df_all <- read.csv(file = file.path("C:/CapProject/25Summer/Model_Result",
                                    paste0(focal_species_under, "_models_result.csv"))
)

df_all_sf <- df_all %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

model_cols <- c("null","LIN","TT","PTT", "M1", "AT", "SQ", "PA")

for (m in model_cols) {
  df_all_sf[[paste0("resid_", m)]] <- df_all_sf$Observed - df_all_sf[[m]]
}

plot_quercus_residual_facet <- function(df_sf, residual_col = "resid_TT") {
  data_model <- df_all_sf
  
  ggplot() +
    ggthemes::theme_few() +
    theme(
      panel.background = element_rect(fill = "gray94", color = NA),
      plot.title = element_text(hjust = 0.5)
    ) +
    
    geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +
    geom_sf(data = nyc_boundary, fill = "white", color = "black") +
    
    geom_sf(data = data_model, aes_string(color = residual_col), size = 1.3, alpha = 0.8) +
    scale_color_gradientn(colors = c("red3", "yellow", "darkgreen"))+
    labs(
      title = paste(focal_species, " - Residuals Map (", gsub("resid_", "", residual_col), ")"),
      color = "Days")+
    annotation_scale(location = "br") +
    facet_wrap(~ Year)
}

res_map <- plot_quercus_residual_facet(df_all_sf, residual_col = paste0("resid_", best_model))

ggsave(
  filename = paste0(focal_species_under, "_res_map.png"),
  plot = res_map,  width = 12,  height = 8,  dpi = 600)


#### Observed Map ####

df_all_sf <- tree_with_pixel %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


plot_quercus_map <- function(df_sf, value_col = "SOS_50") {
  ggplot() +
    ggthemes::theme_few() +
    theme(
      panel.background = element_rect(fill = "gray94", color = NA),
      plot.title = element_text(hjust = 0.5)
    ) +
    
    geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +
    geom_sf(data = nyc_boundary, fill = "white", color = "black") +
    
    geom_sf(data = df_sf, aes_string(color = value_col), size = 1.3, alpha = 0.8) +
    scale_color_gradientn(colors = c("red3", "yellow", "darkgreen")) +
    labs(
      title = paste(focal_species, " - Observed SOS_50"),
      color = "Day of Year"
    ) +
    annotation_scale(location = "br") +
    facet_wrap(~ Year)
}

obs_map <- plot_quercus_map(df_all_sf, value_col = "SOS_50")
print(obs_map)

ggsave(
  filename = paste0(focal_species_under, "_obs_map.png"),
  plot = obs_map,  width = 12,  height = 8,  dpi = 600)



#### All species in same year ####
df_all_species <- bind_rows(df_alba, df_robur, df_rubra, df_velutina)

df_all_species <- df_all_species %>%
  mutate(Residual = Observed - Predicted)

df_all_sf <- df_all_species %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

df_2023_sf <- df_all_sf %>%
  filter(Year == 2023)

# Predicted
pred_2023 <- ggplot() +
  ggthemes::theme_few() +
  theme(
    panel.background = element_rect(fill = "gray94", color = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +
  geom_sf(data = nyc_boundary, fill = "white", color = "black") +
  geom_sf(data = df_2023_sf, aes(color = Predicted), size = 1.3, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("red3", "yellow", "darkgreen"),
    breaks = seq(112, 126, by = 3),  # Predicted values' range
    labels = scales::number_format(accuracy = 1),
    name = "DOY"
  ) +
  facet_wrap(~ Species, ncol = 2) +
  annotation_scale(location = "br")

ggsave(
  filename = "All_Species_pre_2023.png",
  plot = pred_2023,  width = 12,  height = 8,  dpi = 600)



# Residual
resid_2023 <- ggplot() +
  ggthemes::theme_few() +
  theme(
    panel.background = element_rect(fill = "gray94", color = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_sf(data = nyc_boundary_invert, fill = "gray94", color = "gray94") +
  geom_sf(data = nyc_boundary, fill = "white", color = "black") +
  geom_sf(data = df_2023_sf, aes(color = Residual), size = 1.3, alpha = 0.8) +
  scale_color_gradient2(
    low = "red3", mid = "yellow", high = "darkgreen",
    midpoint = 0,
    breaks = seq(-10, 10, by = 5),
    labels = scales::number_format(accuracy = 1),
    name = "Days"
  ) +
  facet_wrap(~ Species, ncol = 2) +
  annotation_scale(location = "br")


ggsave(
  filename = "All_Species_res_2023.png",
  plot = resid_2023,  width = 12,  height = 8,  dpi = 600)
