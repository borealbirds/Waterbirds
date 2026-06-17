# Source dataset DU_SFMN_Wetland_Surveys
# Author: "Melina Houle"
# Date: "March 11, 2026"
## Note on translation:
#     -  
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
library(DBI)
library(odbc)
library(ggplot2)

# config.R store Google Drive ccredential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

organization <- "BU"
dataset_code <- "DU_SFMN"

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
  pid <- "https://drive.google.com/drive/folders/1dG-pFz3GWVnkyceUrbPQInIFgO06rJpf"
  gd.list <- drive_ls(as.character(pid))
  survey <- gd.list %>%
    filter(name =="DU_SFMN_Wetland_Surveys.mdb") %>%
    select("id")
  drive_download(as_id(as.character(survey)), path = file.path(project_dir, "DU_SFMN_Wetland_Surveys.mdb"))
}


con <- dbConnect(odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};", "DBQ=", file.path(project_dir, "DU_SFMN_Wetland_Surveys.mdb"), ";"))
dbListTables(con)

loc <- dbReadTable(con, "EMB Spatial") 

visit <- dbReadTable(con, "tblLandscapeMain") %>%
  dplyr::select(lngLandscapeMain_ID, lngBlock_ID, dteDate, dteTime, lngObserver_ID) %>%
  dplyr::left_join(
    dbReadTable(con, "tblBlock") %>%
      dplyr::select(lngBlock_ID, lngBlock),
    by = "lngBlock_ID"
)

lu_species <- dbReadTable(con, "tblSpeciesList") 

obs <- dbReadTable(con, "tblLandscapeBird") %>% dplyr::select(lngLandscapeMain_ID, lngSpecies_id, intCount, intFirstFiveLastFive, lngDistance_ID, intVisualTime)

survey <- obs %>%
  left_join(visit, by = "lngLandscapeMain_ID") %>% 
  inner_join(loc, by = c("lngBlock" = "Landid"))  %>%
  left_join(lu_species %>% dplyr::select(lngSpecies_ID, strSpeciesCommonName), by = c("lngSpecies_id"="lngSpecies_ID")) %>% 
  dplyr::filter(!is.na(lngSpecies_id)) %>%
  dplyr::filter(!is.na(intCount)) 

# Validate XY _UTM12
xy_sf <- st_as_sf(survey, coords = c("X_coord", "Y_coord"))
xy_sf <- st_set_crs(xy_sf, 3402) %>% st_transform(4326)
canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
bnd <- st_transform(canada, crs= st_crs(4326))
ggplot() +
  geom_sf(data = bnd, fill = NA) +
  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
  theme_minimal()

# test on species
survey_valid <- survey %>%
  dplyr::distinct(strSpeciesCommonName) %>%
  left_join(WT_spTbl, by = c("strSpeciesCommonName"="species_common_name"))

invalid_codes <- survey_valid %>%
  filter(is.na(species_code)) %>%
  distinct(strSpeciesCommonName)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- survey %>%
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
         survey_distance_method	= "0m-50m-100m-INF",
         survey_duration_method	= "0-3-5-10min",
         duration = "10min",
         species_code	= case_when(strSpeciesCommonName == "Unidentified Scaup"  ~ "UNSC",
                                  strSpeciesCommonName == "Western Palm Warbler"  ~ "PAWA",
                                  strSpeciesCommonName == "Black-and-White Warbler" ~ "BAWW",
                                  strSpeciesCommonName == "Gray Jay"  ~ "CAJA",
                                  strSpeciesCommonName == "Slate-colored Junco"  ~ "JUNHYE",
                                  strSpeciesCommonName == "Three-toed Woodpecker"  ~ "ATTW",
                                  strSpeciesCommonName == "Solitary Vireo"  ~ "SOVI",
                                  strSpeciesCommonName == "House Wren"  ~ "NHWR",
                                  strSpeciesCommonName == "Le Conte's Sparrow"  ~ "LCSP",
                                  strSpeciesCommonName == "Unidentified Sandpiper" ~ "UPEE",
                                  strSpeciesCommonName == "Unidentified Dabbler"  ~ "UDAB",
                                  strSpeciesCommonName == "Nelson's Sharp-tailed Sparrow"  ~ "NESP",
                                  strSpeciesCommonName == "Sharp-tailed Sparrow"  ~ "UNSS",
                                  strSpeciesCommonName == "Bald Eagle"  ~ "BAEA",
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


WT.main.report <- rbind(survey_wdur.report, survey_wtime.report, survey_wdist.report, survey_noprot.report)


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
