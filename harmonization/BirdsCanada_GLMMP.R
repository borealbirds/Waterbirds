# Source dataset Birds Canada Great Lakes Marsh Monitoring Program (GLMMP)
# Author: "Melina Houle"
# Date: "February 23, 2026"
## Note on translation:
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
library(naturecounts)

# config.R store Google Drive ccredential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

# Load NatureCount species list and drop all non birds
nc_sp_code <- search_species() %>% dplyr::filter(taxon_group= "BIRDS")

organization <- "BirdsCanada"
dataset_code <- "GLMMP"

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

  pid_glmmp <- "https://drive.google.com/drive/u/0/folders/1W0hCa4kgnuNobNa-xFfkc4afP9ABD22d"
  gd.list <- drive_ls(as.character(pid_glmmp))
  survey_glmmp <- gd.list %>%
    filter(name =="GLMMP.Rdata") %>%
    select("id")
  drive_download(as_id(as.character(survey_glmmp)), path = file.path(project_dir, "GLMMP.Rdata"))
}

load(file.path(project_dir, "GLMMP.Rdata"))
write.csv(glmmp, "./out/BirdsCanada_GLMMP/raw.csv")

survey <- glmmp %>% 
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount)) %>%
  dplyr::filter(!is.na(DecimalLatitude)) %>%
  dplyr::inner_join(sp_Tbl, by = "species_id")  %>%# drop all non-birds
  dplyr::select(GlobalUniqueIdentifier, SamplingEventIdentifier, CollectionCode, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
                CollectorNumber, SurveyAreaIdentifier, TimeObservationsStarted, EffortMeasurement2, ObservationCount, ObservationCount2, 
                ObservationCount3, ObservationCount4, ObservationCount5, ObservationCount6, ObservationCount7, ObservationCount8, ObservationCount9, ObservationCount10,
                ObservationCount11, ScientificName, CommonName, SpeciesCode)

# Validate XY
library(ggplot2)
xy_sf <- st_as_sf(glmmp, coords = c("DecimalLongitude", "DecimalLatitude"))
xy_sf <- st_set_crs(xy_sf, 4326)
canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
bnd <- st_transform(canada, crs= st_crs(4326))
ggplot() +
  geom_sf(data = bnd, fill = NA) +
  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
theme_minimal()

# test on species
glmmp_valid <- glmmp %>%
  dplyr::distinct(SpeciesCode, CommonName) %>%
  left_join(WT_spTbl, by = c("CommonName" ="species_common_name"))

invalid_codes <- glmmp_valid %>%
  filter(is.na(species_code)) %>%
  distinct(CommonName)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- glmmp %>%
  left_join(WT_spTbl %>% select(species_code, scientific_name), by = "scientific_name") %>%
  mutate(organization = "BirdsCanada",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), SurveyAreaIdentifier, sep= ":"),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = DecimalLatitude,
         longitude = DecimalLongitude, 
         start_time = as.POSIXct(as.numeric((TimeObservationsStarted[!is.na(TimeObservationsStarted)])) * 3600, origin = "1970-01-01", tz = "UTC"),
         end_time = as.POSIXct(as.numeric((TimeObservationsEnded[!is.na(TimeObservationsEnded)])) * 3600, origin = "1970-01-01", tz = "UTC"),
         survey_id	= NA,
         survey_time = format(start_time, "%H:%M:%S"),
         survey_date	= paste0(YearCollected, "-", sprintf("%02d", as.integer(MonthCollected)), "-", sprintf("%02d", as.integer(DayCollected)),  " ", survey_time),
         survey_url	= NA,
         observer	= "obsNA",
         protocol_type = "Broadcast Point count",
         survey_distance_method	= "0m-100m-INF/FLY",
         #survey_duration_method	= "0-5-10-15min",
         #duration = round(as.numeric(difftime(end_time, start_time, units = "mins")),0),
         #survey_duration	= case_when(duration <0 ~ "UNKNOWN",
         #                            TRUE ~ paste0(as.character(duration),"min")),
         detection_distance	= case_when(EffortMeasurement2 == "<100m" ~ "0m-100m",
                                        EffortMeasurement2 == "100m"~ "100m-250m",
                                        EffortMeasurement2 == "250m"~ "250m-500m",
                                        EffortMeasurement2 == "500m"~ "500m-1000m",
                                        EffortMeasurement2 == "1 km"~ "1000m-INF",
                                        is.na(EffortMeasurement2) ~ "UNKNOWN"), 
         detection_time	= "UNKNOWN",
         species_code	= case_when(CommonName == "Bald Eagle"  ~ "BAEA",
                                  CommonName == "Eastern Warbling Vireo"  ~ "WAVI",
                                  CommonName == "Northern Yellow Warbler" ~ "YEWA",
                                  CommonName == "Rock Pigeon (Feral Pigeon)"  ~ "ROPI",
                                  CommonName == "moorhen/coot/gallinule sp."  ~ "UNGA",
                                  CommonName == "American Herring Gull"  ~ "HERG",
                                  CommonName == "duck sp."  ~ "UNDU",
                                  CommonName == "gull sp."  ~ "UNGU",
                                  CommonName == "woodpecker sp."  ~ "UNWO",
                                  CommonName == "swallow sp."  ~ "UNSW",
                                  CommonName == "Alder/Willow Flycatcher (Traill's Flycatcher)"  ~ "UNFL",
                                  CommonName == "American Crow (Northwestern)"  ~ "AMCR",
                                  CommonName == "Hudsonian Whimbrel" ~ "WHIM",
                                  CommonName == "new world flycatcher sp." ~ "UNFL",
                                  CommonName == "Eastern Screech-Owl"  ~ "EASO",
                                  CommonName == "swan sp."  ~ "USWN",
                                  CommonName == "Green-winged Teal (American)"  ~ "ANACRE",
                                  CommonName == "Western House-Martin" ~ "COHM",
                                  CommonName == "scoter sp." ~ "USCT",
                                  CommonName == "Spotted/Eastern Towhee (Rufous-sided Towhee)"  ~ "USET",
                                  CommonName == "Bullock's/Baltimore Oriole"  ~ "UNOR",
                                  CommonName == "Scolopacidae sp."  ~ "UNSH",
                                  CommonName == "blackbird sp."  ~ "UNBL",
                                  CommonName == "new world sparrow sp." ~ "UNSP",
                                  CommonName == "bird sp." ~ "UNBI",
                                  CommonName == "Western Flycatcher (Cordilleran)"  ~ "WEFL",
                                  CommonName == "solitary vireo sp." ~ "SOVI",
                                  CommonName == "Redpoll (flammea)" ~ "UNRE",
                                  CommonName == "vireo sp."  ~ "UNVI",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= scientific_name,
         ind_count	= as.numeric(ObservationCount),
         detection_heard	= "f",
         detection_seen	= "t",
         detection_comments = ifelse(english_name %in% c("alcid sp.", "shearwater sp."), paste0("Original species: ", english_name), NA)
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
WT.main.report <- survey %>%
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
save(WT.main.report, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

# Create sub folder in 'toUpload' with the organization name
dr_script<- drive_get("Waterfowl/scripts", shared_drive = "BAM_AvianData")

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
drive_upload(media = script_path, path = as_id(dr_script), name = paste0(organization,"_",dataset_code,".R"), overwrite = TRUE) 
