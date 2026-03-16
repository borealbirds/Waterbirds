# ---
# title: Waterbird models - Data Wrangling
# author: Elly Knight
# created: March 8, 2026
# ---

#NOTES################################

#PURPOSE: This script filters the larger watebird dataset to just point surveys for subsequent steps of analysis.

#This script will require using the locations csv to extract habitat covariates from GEE before running the rest of the code

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"

#3. Load waterbird data object ----
load(file.path(root, "data", "01_WaterBirdData_Wide.Rdata"))

#SUMMARY FOR REPORT #################

#1. Surveys ----
nrow(use)

#2. Families table ---
#This takes a long time! Don't run every time!

# use_fam <- use |> 
#   pivot_longer(LBDO:CITE, values_to="count", names_to="species_code", values_drop_na=TRUE) |> 
#   dplyr::filter(count > 0) |> 
#   left_join(spp |> 
#               dplyr::select(species_code, Family))
# 
# families <- use_fam |> 
#   mutate(survey = case_when(method %in% c("1SPM", "1SPM Audio/Visual hybrid", "1SPT") ~ "ARU",
#                             method %in% c("eBird", "PC", "PC - MMP") ~ "Point count",
#                             method=="Aerial transect" ~ "Aerial transect",
#                             method=="Plot search" ~ "Plot search")) |> 
#   group_by(Family, survey) |> 
#   summarize(detections = n()) |> 
#   ungroup() |> 
#   mutate(percent = round(detections/sum(detections)*100,2)) |> 
#   dplyr::select(-detections) |> 
#   pivot_wider(names_from=survey, values_from = percent, values_fill=0)
#   
# write.csv(families, file.path(root, "results", "FamilySummary.csv"), row.names=FALSE)

#FILTER #############

#1. Remove surveys we don't want ----
dat <- use |> 
  dplyr::filter(!method %in% c("Plot search", "Aerial transect")) |> 
  unique() |> 
  mutate(id = row_number())

#2. Get the locations for covariates ----
locations <- dat |> 
  dplyr::select(location_id, longitude, latitude) |> 
  unique()

#3. Save ----
write.csv(locations, file.path(root, "data", "reporting_mar26", "01_Locations.csv"), row.names=FALSE)

#Use locations to run habitat cov data extraction in GEE and then return to add to this script#

#HABITAT COVS ########

#1. Lookups ----
wc_class <- data.frame(wc_class = c(0:25, 29:32),
                       wetland_og = c("upland",
                                      "lake_freshwater",
                                      "lake_saline",
                                      "reservoir",
                                      "river",
                                      "estuary",
                                      "other",
                                      "stream",
                                      "lacustrine_forest",
                                      "lacustrine_nonforest",
                                      "river_annual_forest",
                                      "river_annual_nonforest",
                                      "river_seasonal_forest",
                                      "river_seasonal_nonforest",
                                      "river_saturated_forest",
                                      "river_saturated_nonforest",
                                      "palustrine_regular_forest",
                                      "palustrine_regular_nonforest",
                                      "palustrine_seasonal_forest",
                                      "palustrine_seasonal_nonforest",
                                      "ephemeral_forest",
                                      "ephemeral_nonforest",
                                      "peatland_arctic_forest",
                                      "peatland_arctic_nonforest",
                                      "peatland_temperate_forest",
                                      "peatland_tempoerate_nonforest",
                                      "saltmasch",
                                      "delta",
                                      "other_coastal",
                                      "saltpan")) |> 
  mutate(wetland = case_when(str_detect(wetland_og, "river") ~ "river",
                             str_detect(wetland_og, "palustrine") ~ "palustrine",
                             str_detect(wetland_og, "lacustrine") ~ "lacustrine",
                             str_detect(wetland_og, "peatland") ~ "peatland",
                             str_detect(wetland_og, "ephemeral") ~ "ephemeral",
                             str_detect(wetland_og, "lake") ~ "lake",
                             wetland_og %in% c("saltmarsh", "delta", "other_coastal", "saltpan") ~ "coastal",
                             !is.na(wetland_og) ~ wetland_og),
         forested = ifelse(str_detect(wetland_og, "_forest"), "forested", "nonforest")) |> 
  mutate(wetland = factor(wetland),
         wetland = relevel(wetland, "upland"),
         forested = factor(forested, levels = c("nonforest", "forested")))

#2. Get the water data ----
water_2 <- read.csv(file.path(root, "data", "covariates", "WaterbirdLocations_water_2km.csv"))
water_20 <- read.csv(file.path(root, "data", "covariates", "WaterbirdLocations_water_2km.csv"))

#3. Get the wetland class data ----
wc <- read.csv(file.path(root, "data", "covariates", "WaterbirdLocations_wetlandclass_point.csv"))

#5. Wrangle covariates ----
cov <- dat |> 
  dplyr::select(id, location_id, longitude, latitude) |> 
  left_join(wc |> 
               rename(wc_class = mean) |> 
               dplyr::select(location_id, wc_class)) |> 
  left_join(water_2 |> 
               rename(uplandp_2km = b1,
                      wetp_2km = transition,
                      occurrence_2km = occurrence,
                      recurrence_2km = recurrence,
                      seasonality_2km = seasonality) |> 
               dplyr::select(location_id, uplandp_2km, occurrence_2km, recurrence_2km, seasonality_2km, wetp_2km)) |> 
  left_join(water_20 |> 
              rename(uplandp_20km = b1,
                     wetp_20km = transition,
                     occurrence_20km = occurrence,
                     recurrence_20km = recurrence,
                     seasonality_20km = seasonality) |> 
              dplyr::select(location_id, uplandp_20km, occurrence_20km, recurrence_20km, seasonality_20km, wetp_20km)) |> 
  left_join(wc_class) |> 
  mutate(occurrence_2km = ifelse(is.na(occurrence_2km), 0, 1),
         recurrence_2km = ifelse(is.na(recurrence_2km), 0, 1),
         seasonality_2km = ifelse(is.na(seasonality_2km), 0, 1),
         occurrence_20km = ifelse(is.na(occurrence_20km), 0, 1),
         recurrence_20km = ifelse(is.na(recurrence_20km), 0, 1),
         seasonality_20km = ifelse(is.na(seasonality_20km), 0, 1))

#6. Save ----
save(dat, cov, spp, file=file.path(root, "data", "reporting_mar26", "01_WaterBirdData_Habitat.Rdata"))
