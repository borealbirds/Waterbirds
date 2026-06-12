# ---
# title: Waterbird models - Data - Wrangle Dataset ----
# author: Elly Knight
# created: March 8, 2026
# ---

#NOTES################################

#PURPOSE: This script makes a figure of the waterbird dataset for the 2025/26 BAM annual report.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"

#3. Get data ----
load(file.path(root, "data", "01_WaterBirdData_Wide.Rdata"))

#MAKE A FIGURE #############

#1. Get Canada boundary ----
can <- read_sf(file.path("G:/Shared drives/BAM_NationalModels5", "Regions", "CAN_adm", "CAN_adm1.shp")) |>
  st_transform(crs=4326)

#1. Aggregate rows to location ----
locs <- use |> 
  mutate(lonr = round(longitude, 2),
         latr = round(latitude, 2),
         year = year(date_time),
         method = case_when(method %in% c("1SPM", "1SPM Audio/Visual hybrid", "1SPT") ~ "ARU",
                            method=="PC - MMP" ~ "Broadcast survey",
                            method=="Aerial transect" ~ "Aerial survey",
                            method=="PC" ~ "Point count",
                            !is.na(method) ~ method)) |> 
  group_by(method, lonr, latr) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  st_as_sf(coords=c("lonr", "latr"), crs=4326, remove=FALSE) |> 
  st_intersection(can) |> 
  mutate(surveys = ifelse(n > 50, 50, n))

#3. Plot ----
method.plot <- ggplot() +
  geom_sf(data=can, fill=NA, colour="black") +
  geom_point(data=locs, aes(x=lonr, y=latr, colour=surveys), alpha = 0.5, size=0.2) +
  scale_colour_viridis_c(name="Number of surveys") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=10)) +
  facet_wrap(~method)

ggsave(method.plot, file=file.path(root, "figures", "Datasetdistribution.jpeg"), width = 12, height=8)
