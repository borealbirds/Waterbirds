# Source dataset Birds Canada BC Marsh Monitoring Program (BCMMP)
# Author: "Melina Houle"
# Date: "February 23, 2026"
## Note on translation:
#     - join NC_species table to join on data to drop all non-bird species
#     - Drop ObservationCount12 = "1" (Before/After)
#     - Only 400 /16 000 obs have distance band. Distance band is used as TRUE/FALSE. If TotalCount =12 and observationCCount9 (<50m) = 1 and 
#       ObservationCount10 (50-100m) =1, there is no way to retrieve how many bird were seen in each distance band. I marked them as UNKNOWN. If only one distance band is TRUE, 
#       I used the band. 
#     - TotalCount is splitted by duration band most of the time.  
#     - When observationCount3 (0-3min), ObservationCount4 (3-5min) and ObservationCount5 (5-10min) are all 0 but ObservationCount2 (first minute detected) is >0,
#       I used that to attribute duration band. 
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
library(reshape2)

# config.R store Google Drive ccredential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

# NatureCount species
#install.packages("naturecounts", 
#                 repos = c(birdscanada = 'https://birdscanada.r-universe.dev',
#                           CRAN = 'https://cloud.r-project.org'))
library(naturecounts)
sp_Tbl <- search_species(show = "all") %>%
  filter(class== "Aves")

organization <- "BirdsCanada"
dataset_code <- "BCMMP"

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

  pid <- "https://drive.google.com/drive/u/0/folders/1W0hCa4kgnuNobNa-xFfkc4afP9ABD22d"
  gd.list <- drive_ls(as.character(pid))
  survey <- gd.list %>%
    filter(name =="BCMMP.Rdata") %>%
    select("id")
  drive_download(as_id(as.character(survey)), path = file.path(project_dir, "BCMMP.Rdata"))
  
}

load(file.path(project_dir, "BCMMP.Rdata"))

bcmmp_reduced <- bcmmp %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount) & ObservationCount12 != "1") %>%
  dplyr::filter(!is.na(DecimalLatitude)) %>%
  dplyr::inner_join(sp_Tbl, by = "species_id")  %>%# drop all non-birds
  dplyr::select(GlobalUniqueIdentifier, SamplingEventIdentifier, CollectionCode, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
                CollectorNumber, SurveyAreaIdentifier, TimeObservationsStarted, EffortMeasurement2, ObservationCount, ObservationCount2, 
       ObservationCount3, ObservationCount4, ObservationCount5, ObservationCount6, ObservationCount7, ObservationCount8, ObservationCount9, ObservationCount10,
       ObservationCount11, ScientificName, CommonName, SpeciesCode)

# Validate XY
#xy_sf <- st_as_sf(accws, coords = c("DecimalLongitude", "DecimalLatitude"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()

# test on species
bcmmp_valid <- bcmmp_reduced %>%
  dplyr::distinct(SpeciesCode, ScientificName, CommonName) %>%
  left_join(WT_spTbl, by = c("ScientificName"="scientific_name"))

invalid_codes <- bcmmp_valid %>%
  filter(is.na(species_code)) %>%
  distinct(CommonName)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- bcmmp_reduced %>%
  left_join(WT_spTbl %>% select(species_code, scientific_name), by = c("ScientificName"="scientific_name")) %>%
  mutate(organization = "BirdsCanada",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), SurveyAreaIdentifier, sep= ":"),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = DecimalLatitude,
         longitude = DecimalLongitude, 
         start_time = as.POSIXct(as.numeric((TimeObservationsStarted[!is.na(TimeObservationsStarted)])) * 3600, origin = "1970-01-01", tz = "UTC"),
         survey_id	= NA,
         survey_time = format(start_time, "%H:%M:%S"),
         survey_date	= paste0(YearCollected, "-", sprintf("%02d", as.integer(MonthCollected)), "-", sprintf("%02d", as.integer(DayCollected)),  " ", survey_time),
         survey_url	= NA,
         observer	= paste0("obs", as.character(CollectorNumber)),
         protocol_type = "Stationary Count Broadcast",
         survey_distance_method	= "0m-50m-100m-INF/FLY",
         survey_duration_method	= "0-3-5-10min",
         duration = "10min",
         species_code	= case_when(CommonName == "Cooper's Hawk"  ~ "COHA",
                                  CommonName == "Hudsonian Whimbrel"  ~ "WHIM",
                                  CommonName == "Rock Pigeon (Feral Pigeon)" ~ "ROPI",
                                  CommonName == "Hairy Woodpecker"  ~ "HAWO",
                                  CommonName == "Northern Yellow Warbler"  ~ "YEWA",
                                  CommonName == "Evening Grosbeak"  ~ "EVGR",
                                  CommonName == "Western Warbling Vireo"  ~ "WAVI",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)]
         )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

# Four scenarios ( n = 16625)
## 1 - Obs has duration band
## 2 - Obs has no duration band but 1st min detected
## 3 - obs has no duration band but distance band. 
## 4 - obs has no info on distance or duration
# 1. Obs has duration band (n = 14177)
survey_wdur <- survey %>%
  dplyr::filter(dplyr::if_any(c(ObservationCount3, ObservationCount4, ObservationCount5), ~ . != "0")) %>%
  reshape2::melt(measure.vars = c("ObservationCount3","ObservationCount4","ObservationCount5"), value.name = "abundance")  %>%
  dplyr::filter(abundance != 0) %>%
  mutate(survey_duration	= case_when(variable == "ObservationCount3" ~ "0-3min",
                                     variable == "ObservationCount4" ~ "3-5min",
                                     variable == "ObservationCount5" ~ "5-10min"),
         detection_distance	= case_when(ObservationCount9 == "1" & ObservationCount10 =="0" & ObservationCount11 =="0" ~ "0m-50m",
                                        ObservationCount9 == "0" & ObservationCount10 =="1" & ObservationCount11 =="0" ~ "50m-100m",
                                        ObservationCount9 == "0" & ObservationCount10 =="0" & ObservationCount11 =="1" ~ "100m-INF",
                                        ObservationCount6 == "1" ~ "100m-INF/FLY",
                                        TRUE ~ "UNKNOWN"),
         detection_time	= survey_time,
         ind_count	= as.numeric(abundance),
         detection_heard	= case_when(ObservationCount7 =="1"~ "t",
                                     TRUE ~ "dnc"),
         detection_seen	= case_when(ObservationCount8 =="1"~ "t",
                                    TRUE ~ "dnc"),
         detection_comments = NA)


# 2. Obs has no duration band but 1st min detected (n = 812)
#    - 4 obs have multiple distance band. They all have a count of 1. I assigned the closest distance band to the obs. 
#    - Some obs have 0 for aural and visual. The columncs can't be interpreted as true FALSE. 
survey_wtime <- survey %>%
  dplyr::filter(dplyr::if_all(c(ObservationCount3, ObservationCount4, ObservationCount5), ~ . == "0") & ObservationCount2 != "0") %>%
  mutate(survey_duration	= case_when(ObservationCount2 %in% c("1", "2", "3") ~ "0-3min",
                                     ObservationCount2 %in% c("4", "5") ~ "3-5min",
                                     ObservationCount2 %in% c("6", "7", "8", "9","10") ~ "5-10min"),
         detection_distance	= case_when(ObservationCount10 == "1" ~ "50m-100m",
                                        ObservationCount9 == "1" ~ "0m-50m",
                                        ObservationCount11 == "1" ~ "100m-INF"),
         detection_time	= survey_time,
         ind_count	= as.numeric(ObservationCount),
         detection_heard	= case_when(ObservationCount7 =="1"~ "t",
                                     TRUE ~ "dnc"),
         detection_seen	= case_when(ObservationCount8 =="1"~ "t",
                                    TRUE ~ "dnc"),
         detection_comments = NA)


# 3 - obs has no duration band but a known distance band. (n = 16). No need to pivot. Distance band always unique.
#     all obs are either classified as aural (ObservationCOunt7) or visual (ObservationCount8)
survey_wdist <- survey %>%
  dplyr::filter(dplyr::if_all(c(ObservationCount3, ObservationCount4, ObservationCount5), ~ . == "0") & ObservationCount2 == "0"
                & dplyr::if_any(c(ObservationCount9, ObservationCount10, ObservationCount11), ~ . != "0")) %>%
  mutate(survey_duration	= "UNKNOWN",
         detection_distance	= case_when(ObservationCount9 == "1" ~ "0m-50m",
                                        ObservationCount10 == "1" ~ "50m-100m",
                                        ObservationCount11 == "1" ~ "100m-INF"),
         detection_time	= survey_time,
         ind_count	= as.numeric(ObservationCount),
         detection_heard	= case_when(ObservationCount7 =="1"~ "t",
                                     TRUE ~ "f"),
         detection_seen	= case_when(ObservationCount8 =="1"~ "t",
                                    TRUE ~ "f"),
         detection_comments = NA)

# 4 - obs has no duration and distance band. (n = 1620). ObservationCount6 means outside/flyover. 
survey_noprot <- survey %>%
  dplyr::filter(dplyr::if_all(c(ObservationCount3, ObservationCount4, ObservationCount5), ~ . == "0") & ObservationCount2 == "0"
                & dplyr::if_all(c(ObservationCount9, ObservationCount10, ObservationCount11), ~ . == "0")) %>%
  mutate(survey_duration	= "UNKNOWN",
         detection_distance	= case_when(ObservationCount6 == "1" ~ "100m-INF/FLY",
                                        TRUE ~ "UNKNOWN"),
         detection_time	= survey_time,
         ind_count	= as.numeric(ObservationCount),
         detection_heard	= case_when(ObservationCount7 =="1"~ "t",
                                     TRUE ~ "dnc"),
         detection_seen	= case_when(ObservationCount8 =="1"~ "t",
                                    TRUE ~ "dnc"),
         detection_comments = NA)

        
        

#---
survey_wdur.report <- survey_wdur %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, survey_duration, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")

survey_wtime.report <- survey_wtime %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, survey_duration, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")

survey_wdist.report <- survey_wdist %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, survey_duration, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")

survey_noprot.report <- survey_noprot %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, survey_duration, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")


BCMMP <- rbind(survey_wdur.report, survey_wtime.report, survey_wdist.report, survey_noprot.report)


#--------------------------------------------------------------
#
#       EXPORT
#
#--------------------------------------------------------------

# Create sub folder in 'toUpload' with the organization name
dr<- drive_get("waterfowl_Rdata/", shared_drive = "BAM_AvianData")

# Save .Rdata
save(BCMMP, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 


## Save WT upload location 
location <- BCMMP %>%
  dplyr::select(organization, location, longitude, latitude, location_buffer_m) %>%
  dplyr::rename(buffer_m = location_buffer_m) %>%
  mutate(location_visibility = "Visible",
         true_coordinates = TRUE,
         location_comments = NA,
         internal_wildtrax_id = NA)

write.csv(location, file= file.path(out_dir, paste0(organization,"_", dataset_code, "_location.csv")), row.names = FALSE, na = "")
location_out <- file.path(out_dir, paste0(organization,"_", dataset_code,"_location.csv"))
drive_upload(media = location_out, path = as_id(dr), name = paste0(organization,"_", dataset_code,"_location.csv"), overwrite = TRUE) 
  
## Save WT upload survey
survey <- BCMMP %>%
  dplyr::select(location, survey_date, survey_duration_method, survey_distance_method, observer, species_code, detection_distance,
                survey_duration, individual_count, detection_seen, detection_heard, detection_comments) %>%
  dplyr::rename(durationMethod = survey_duration_method,
                distanceMethod = survey_distance_method,
                species = species_code, 
                distanceband = detection_distance,
                durationinterval= survey_duration, 
                abundance = individual_count,
                isHeard = detection_heard,
                isSeen = detection_seen,
                comments = detection_comments)

write.csv(survey, file= file.path(out_dir, paste0(organization,"_", dataset_code, "_survey.csv")), row.names = FALSE, na = "")
survey_out <- file.path(out_dir, paste0(organization,"_", dataset_code,"_survey.csv"))
drive_upload(media = survey_out, path = as_id(dr), name = paste0(organization,"_", dataset_code,"_survey.csv"), overwrite = TRUE) 

## Save script
script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
