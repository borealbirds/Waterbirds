# Source dataset USWFS Waterfowl Breeding Population and Habitat Survey (WBPHS)
# Author: "Melina Houle"
# Date: "February 2, 2026"
## Note on translation:
#      -species dictionary is found in WBPHS_Geolocated_Counts.html 
#----------------------------------------------
#update.packages()
library(dplyr) # mutate, %>%
library(readr) #read_delim
library(stringr) #str_replace_all
library(googledrive) #drive_get, drive_mkdir, drive_ls, drive_upload
library(terra)
library(sf)
library(hms)
library(readxl)

# config.R store Google Drive ccredential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

#species
WT_spTbl <- read.csv(file.path("./lookupTables/species_codes.csv"),fileEncoding = "UTF-8-BOM") %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL"))

organization <- "USFWS"
dataset_code <- "WBPHS"

out_dir <- file.path("./out", paste(organization, dataset_code, sep= "_"))   # where output dataframe will be exported
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}
project_dir <- file.path("./project", paste(organization, dataset_code, sep= "_"))
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}
 
#--------------------------------------------------------------
#
#       DOWNLOAD FILE FROM DRIVE 
#
#--------------------------------------------------------------
projFold <- "WBPHS_Geolocated_Counts"
if (length(list.files(project_dir)) ==0) {
  url <- "https://iris.fws.gov/APPS/ServCat/DownloadFile/277351"

  download.file(url, file.path(project_dir, paste0(projFold, ".zip")), mode = "wb")
  unzip(file.path(project_dir, paste0(projFold, ".zip")), exdir = project_dir)
}

wbphs <- read.csv(file.path(project_dir, projFold, "wbphs_geolocated_counts_forDistribution.csv"), as.is = T, na.strings = c("NULL"), fileEncoding="UTF-8-BOM") %>%
  dplyr::filter(!is.na(latitude)) %>%
  dplyr::filter(survey_species != "POND")

#xy_sf <- st_as_sf(wbphs, coords = c("longitude", "latitude"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#usa <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp")
#bnd <- st_transform(usa, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()


# test on species
wbphs_valid <- wbphs %>%
  distinct(survey_species) %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name),
            by = c("survey_species" = "species_code"))
invalid_codes <- wbphs_valid %>%
  filter(is.na(species_common_name)) %>%
  distinct(survey_species)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- wbphs %>%
  mutate(organization = "USFWS",
         project= dataset_code,
         project_id	= NA,
         location	= paste0(dataset_code, ":s", stratum, "t",transect, "seg", segment),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = latitude,
         longitude	= longitude,
         survey_id	= NA,
         survey_time = "00:00:01",
         survey_date	= paste0(survey_year, "-", sprintf("%02d", as.integer(survey_month)), "-", sprintf("%02d", as.integer(survey_day)), " ", survey_time),
         survey_url	= NA,
         observer	= "obsNA",
         protocol_type = "Segment Transect",
         survey_distance_method	= "UNKNOWN",
         survey_duration_method	= "UNKNOWN",
         survey_duration	= "UNKNOWN",
         detection_distance	= "UNKNOWN",
         detection_time	= "UNKNOWN",
         species_code	= case_when(survey_species  == "SCAU"  ~ "UNSC",
                                  survey_species  == "AGWT"  ~ "GWTE",
                                  survey_species  == "RNDU"  ~ "RNDU",
                                  survey_species  == "GOLD"  ~ "UGOL",
                                  survey_species  == "NOPI" ~ "NOPI",
                                  survey_species  == "NSHO"  ~ "NSHO",
                                  survey_species  == "SWAN"  ~ "USWN",
                                  survey_species  == "SCOT"  ~ "USCT",
                                  survey_species  == "CAGO"  ~ "CANG",
                                  survey_species  == "MERG" ~ "UNME",
                                  survey_species  == "REDH"  ~ "REDH",
                                  survey_species  == "EIDE"  ~ "UDIV",
                                  survey_species  == "RUDU"  ~ "RUDU",
                                  TRUE ~ survey_species),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= count,
         detection_heard	= "dnc",
         detection_seen	= "dnc",
         detection_comments = social_group ) 

#---
WBPHS <- survey %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, survey_duration, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")

# test on duplicates
#dup <- survey %>%
#    dplyr::count(across(everything())) %>%
#    dplyr::filter(n > 1)

#  View(dup)

#--------------------------------------------------------------
#
#       EXPORT
#
#--------------------------------------------------------------

# Create sub folder in 'toUpload' with the organization name
dr<- drive_get("waterfowl_Rdata/", shared_drive = "BAM_AvianData")

# Save
save(WBPHS, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
