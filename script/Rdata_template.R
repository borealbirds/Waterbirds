# Source dataset 
# Author: 
# Date: 
## Note on translation:
#         - obs that were dropped
#         - specification on dataset

#----------------------------------------------
#update.packages()
library(dplyr) # mutate, %>%
library(readr) #read_csv
library(stringr) #str_replace_all
library(googledrive) #drive_get, drive_mkdir, drive_ls, drive_upload
library(googlesheets4)
library(terra)
library(sf)
library(hms)
library(httr)
library(readxl)

# config.R store Google Drive credential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

#species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv"),locale = locale(encoding = "UTF-8")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL"))

WT_durMethTbl <- read.csv(file.path("./lookupTables/duration_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distMethTbl <- read.csv(file.path("./lookupTables/distance_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_durBandTbl <- read.csv(file.path("./lookupTables/duration_interval_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distBandTbl <- read.csv(file.path("./lookupTables/distance_band_codes.csv"), fileEncoding="UTF-8-BOM")

#Observer
obs_url <- "https://docs.google.com/spreadsheets/d/1gsm4LSwU31vJQIh5Ahpy70dhYvPr9ftHqaW1gw75IeU"
observer_Tbl <-  read_sheet(obs_url, sheet = "master_observer.csv")

organization <- "***Organization***"
dataset_code <- "***ProjectCode***"
lu <- "./lookupTables"

out_dir <- file.path("./out", dataset_code)    # where output dataframe will be exported
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}
project_dir <- file.path(wd, "project", dataset_code)
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}
 
#--------------------------------------------------------------
#
#       DOWNLOAD FILE FROM DRIVE 
#
#--------------------------------------------------------------
if (length(list.files(project_dir)) ==0) {
  pid <- "https://drive.google.com/drive/folders/1psX3W0yoeTY2K-UO6TXEWSZfU-speC0f"
  #Download from GoogleDrive
  gd.list <- drive_ls(as.character(pid))
  
  survey <- gd.list %>%
    filter(id =="1JIWayZtJBplltH2kLA9WLRPTeC5brjOG") %>%
    select(id)
  drive_download(as_id(as.character(survey)), path = file.path(project_dir, "filename"))
  
}
  
survey_tbl <- read_excel(file.path(project_dir, "filename.xlsx"), sheet = "Sheet1")

 ## Validate mapping to make sure the point fall where it is expected
 # library(ggplot2)
 # df_poly <- st_as_sf(survey_tbl, coords = c("Longitude", "Latitude"), crs = 4326) 
 # canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
 # bnd <- st_transform(canada, crs= st_crs(4326))
 # ggplot() +
 #   geom_sf(data = bnd, fill = NA) +
 #   geom_sf(data = df_poly, fill = NA, colour = "red", linewidth = 6) +
 #   theme_minimal()

#################################################
## HARMONIZE
#################################################
survey <-  survey_tbl %>%
  dplyr::filter(Latitude >0) %>%
  mutate(organization = "***Organization***",
       project= dataset_code,
       project_id	= NA,
       location	= paste0(dataset_code, ":", Polygon, ":", FID_),
       location_id	= NA,
       location_buffer_m	= NA,
       latitude = Latitude,
       longitude	= Longitude,
       survey_id	= NA,
       survey_date	= as.Date(paste(Year, sprintf("%02d", Month), sprintf("%02d", Day), sep = "-")),
       observer	= "***obsNA***",
       survey_distance_method	= "***AreaSearch_NA***",
       survey_duration_method = "***UNKNOWN***",
       survey_duration	= "***UNKNOWN***",
       detection_distance	= paste0(as.character(round(dist,0)), "m"),
       detection_time	= as.character(format(as.POSIXct(Time, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S")), 
       survey_time = "***UNKNOWN***",
       surveyDateTime = paste(survey_date, survey_time),
       survey_url	= NA,
       ind_count	= as.integer(Count),
       detection_heard	= "***No***",
       detection_seen	= "***Yes***"
       )
  
  ################################
  #### Update master_observer ####
  ################################
  unique_observers <- tibble(
      observer_name = "***obsNA***",
      observer_id = "***obsNA***"
    )
  
  # Create the append_obs data frame
  append_obs <- unique_observers %>%
    select(observer_id, observer_name) %>%
    mutate(
      organization = "***Organization***",
      project = dataset_code
    ) %>%
    select(organization, project, observer_id, observer_name)
  
  # Identify rows in append_obs that are not in observer_Tbl
  new_rows <- anti_join(append_obs, observer_Tbl, 
                        by = c("organization", "project", "observer_id", "observer_name"))
  
  # Combine new rows with the existing observer_Tbl
  if (nrow(new_rows) > 0) {
    sheet_append(obs_url, new_rows)
  }
  
  ###########################################
  ## SPECIES
  ###########################################
  
  #list the ones that won't pass
  unique(survey[!survey$Paul_Spp %in% WT_spTbl$species_code, ]$Paul_Spp)
  
  # Example of manual fix. 
  survey <- survey %>% 
    mutate(species_code = case_when(Paul_Spp  == "CAGO" ~ "CANG",
                               Paul_Spp  == "Eider" ~ "UDIV",
                               Paul_Spp  == "CRANE" ~ "UNBI", 
                               Paul_Spp  == "DUCKS" ~ "UNDU", 
                               Paul_Spp  == "Gull" ~ "UNGU", 
                               Paul_Spp  == "JAEG" ~ "UJAE", 
                               Paul_Spp  == "LOON" ~ "UNLO", 
                               Paul_Spp  == "RLHAWK" ~ "RLHA", 
                               Paul_Spp  == "ROSS" ~ "ROGO", 
                               Paul_Spp  == "TERN" ~ "UNTE", 
                               TRUE ~ WT_spTbl$species_code[match(Paul_Spp, WT_spTbl$species_code)]),
           species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
           species_scientific_name = WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
           detection_comments = case_when(Paul_Spp  == "Eider" ~ "Unidentified eider",
                                          Paul_Spp  == "CRANE" ~ "Unidentified crane",
                                          TRUE ~ NA)
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
  #dup <- WT.main.report %>%
  #  dplyr::count(across(everything())) %>%
  #  dplyr::filter(n > 1)

  #View(dup)

#--------------------------------------------------------------
#
#       EXPORT
#
#--------------------------------------------------------------

# Create sub folder in 'toUpload' with the organization name
dr<- drive_get("waterfowl_Rdata/", shared_drive = "BAM_AvianData")

# Save
save(WT.main.report, file = file.path(out_dir, paste0(dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization, "_", dataset_code,".RData"), overwrite = TRUE) 

dr<- drive_get("Waterfowl/scripts", shared_drive = "BAM_AvianData")

script_path <- file.path(wd, "script", paste0(organization, "_", dataset_code,".R"))
drive_upload(media = script_path, path = as_id(dr), name = paste0(organization, "_", dataset_code,".R"), overwrite = TRUE) 
