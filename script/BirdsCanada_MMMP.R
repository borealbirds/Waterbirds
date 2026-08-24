# Source dataset Birds Canada Maritimes Marsh Monitoring Program (MMMP)
# Author: "Melina Houle"
# Date: "August 24, 2026"
## Note on translation:
#    -ObservationCount ="Total Count"
#    -ObservationCount2 ="0-5 min"
#    -ObservationCount3 ="5-10 min"
#    -ObservationCount4 ="10-15 min"
#    -ObservationCount5 ="Target species (5-15 min)"
#    -ObservationCount6 ="Outside/Flythrough"
#    - Drop 20 obs with unknown SpeciesCode
#  Based on the description of ObservationCount columns, we used totalCount for abundance and try to split by time interval 0-5, 5-10 or 10-15 when known. 
#  We log the flyover but 
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
library(tibblify)
library(reshape2)

install.packages("naturecounts", 
                 repos = c(birdscanada = 'https://birdscanada.r-universe.dev',
                           CRAN = 'https://cloud.r-project.org'))
library(naturecounts)

# config.R store Google Drive credential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

# NatudeCount species
sp_Tbl <- search_species() %>% filter(taxon_group == "BIRDS")

organization <- "BirdsCanada"
dataset_code <- "MMMP"

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
if (length(list.files(project_dir)) ==0) {
  mmp <- nc_data_dl(collections = "MMMP", fields_set = "extended", request_id = r_id, username = natureCounts_username, info = "download Atlas")
  save(mmp, file = file.path(project_dir, "mmpbirds.RData"))
}else{
  load(file.path(project_dir, "mmpbirds.RData"))
}

# test on species
mmp_invalid <- mmp %>%
  dplyr::distinct(species_id, CommonName, SpeciesCode) %>%
  filter(is.na(species_id)) %>%
  distinct(SpeciesCode)

mmp_invalid
mmp_remove <- c("BLBW?", "RUGR?", "GCFL?", "COTE?", "NOWA?", "RBNU?", "ATSP?", "GBHE?", "TRES?", "BLTRBL", "COYT","PHBE", "AM", "NORA", "COYH",
                 "BLCH", "REST", "WHSP", "REVE", "AMPH")
survey <- mmp %>% 
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount)) %>%
  dplyr::filter(!is.na(DecimalLatitude)) %>%
  dplyr::filter(!SpeciesCode %in% mmp_remove) %>% 
  dplyr::select(GlobalUniqueIdentifier, SamplingEventIdentifier, CollectionCode, latitude, longitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
                CollectorNumber, SurveyAreaIdentifier, TimeObservationsStarted, EffortMeasurement1, ObservationCount, ObservationCount2, 
                ObservationCount3, ObservationCount4, ObservationCount5, ObservationCount6, ScientificName, CommonName, SpeciesCode, species_id)

# Validate XY
library(ggplot2)
xy_sf <- st_as_sf(survey, coords = c("longitude", "latitude"))
xy_sf <- st_set_crs(xy_sf, 4326)
canada <- st_read("./lookupTables/GIS_layers/CanadaLAEA.shp")
bnd <- st_transform(canada, crs= st_crs(4326))
ggplot() +
  geom_sf(data = bnd, fill = NA) +
  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
theme_minimal()

# test on species with WildTRAX
mmp_species <- survey %>%
  dplyr::distinct(SpeciesCode, CommonName ) %>%
  left_join(WT_spTbl, c("CommonName" = "species_common_name"))

invalid_codes <- mmp_species %>%
  filter(is.na( species_code)) %>%
  distinct()

invalid_codes
#################################################
## HARMONIZE
#################################################
survey <- survey %>%
  mutate(ObservationDescription7 = "unknown band",
         ObservationCount7 = as.character(as.numeric(ObservationCount) -rowSums(across(c(ObservationCount2, ObservationCount3,
                                                                                         ObservationCount4),
                                                                                       as.numeric), na.rm = TRUE)
                                          ))
  
survey_expanded <- melt(survey, measure.vars = c("ObservationCount2","ObservationCount3","ObservationCount4","ObservationCount6", "ObservationCount7"), value.name = "abundance")
survey_expanded$abundance <- as.numeric(survey_expanded$abundance)
survey_expanded <- subset(survey_expanded, survey_expanded$abundance>0 )

survey.flat <- survey_expanded %>%
  rename(species_code = SpeciesCode) %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name, scientific_name), by = "species_code") %>%
  mutate(organization = "BirdsCanada",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), SurveyAreaIdentifier, sep= ":"),
         location_id	= NA,
         location_buffer_m = NA,
         survey_id	= NA,
         survey_time = ifelse(is.na(TimeCollected), "00:00:01", TimeCollected),
         survey_date	= paste0(YearCollected, "-", sprintf("%02d", as.integer(MonthCollected)), "-", sprintf("%02d", as.integer(DayCollected)),  " ", survey_time),
         survey_url	= NA,
         observer	= "obsNA",
         protocol_type = "Broadcast Point count",
         survey_distance_method	= "0m-100m",
         survey_duration_method	= "0-5-10-15min",
         detection_distance	= "0m-100m", 
         detection_time	= case_when(variable == "ObservationCount2" ~ "0-5min",
                                    variable == "ObservationCount3" ~ "5-10min",
                                    variable == "ObservationCount4" ~ "10-15min",
                                    variable == "ObservationCount6" ~ "flyover",
                                    variable == "ObservationCount7" ~ "UNKNOWN"),
         species_code	= case_when(species_code == "NOWT" ~ "NOWA",
                                  species_code == "SWALLOW SP" ~ "UNSW",
                                  species_code == "PEEP SP"  ~ "UPEE",
                                  species_code == "WARBLER SP"  ~ "UNWA",
                                  species_code == "SANDPIPER SP"  ~ "UPEE",
                                  species_code == "HERG"  ~ "AHGU",
                                  species_code == "BDOW"  ~ "BAOW",
                                  species_code == "GOOT"  ~ "UNGA",
                                  TRUE ~ species_code),
         species_scientific_name	= scientific_name,
         ind_count	= abundance,
         detection_heard	= "dnc",
         detection_seen	= "dnc",
         detection_comments = case_when(species_code == "UNGA" ~ "moorhen/coot/gallinule sp.",
                                        TRUE ~ "Broadcast 0-5min, silent 5-15min")
         )

# check
print(unique(survey.flat$species_code[!(survey.flat$species_code %in% WT_spTbl$species_code)]))

#---
MMP <- survey.flat %>%
  dplyr::filter(!detection_time == "flyover") %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")

#--------------------------------------------------------------
#
#       EXPORT
#
#--------------------------------------------------------------

# Create sub folder in 'toUpload' with the organization name
dr<- drive_get("waterfowl_Rdata/", shared_drive = "BAM_AvianData")

# Save .Rdata
save(MMP, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 
