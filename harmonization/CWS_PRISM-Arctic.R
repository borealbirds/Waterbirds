# Source dataset CWS PRISM-Arctic (Tiers 1)
# Author: "Melina Houle"
# Date: "May 11, 2026"
## Note on translation:
#         - Delete all non birds obs
#         - Delete no bird obs
#         - Wasn't able to download kmz with drive_download. I had to use GET. 
#         - Kept outside plot and used Sighting_Type as detection_distance. 
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

# config.R store Google Drive ccredential and path to working directory
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

organization <- "CWS"
dataset_code <- "PRISM-Arctic"
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
  pid <- "https://drive.google.com/drive/folders/1UJILNvsZfgE9HNutBH3B3R2uARJMg5NB"
  #Download from GoogleDrive
  gd.list <- drive_ls(as.character(pid))
  
  survey <- gd.list %>%
    filter(id =="1TSvFr4OFFRMaZNSTrkZQgXo0jVTrNywH") %>%
    select(id)
  drive_download(as_id(as.character(survey)), path = file.path(project_dir, "Arctic_PRISM_Surveys_Tiers1_Birds_Data_Donnees.csv"))
  
  metadata <- gd.list %>%
    filter(id =="1uvsW50MAmL4YoLpR6_6UFfy8KsPZr7h2") %>%
    select(id)
  drive_download(as_id(as.character(metadata)), path = file.path(project_dir, "Arctic_PRISM_Surveys_Tier1_Metadata.txt"))
  
  species_lu <- gd.list %>%
    filter(id =="1BPnIy3otp7whZaTonus2uwVBqHzv1Mph") %>%
    select(id)
  drive_download(as_id(as.character(species_lu)), path = file.path(project_dir, "Arctic_PRISM_Species_Code.csv"))
  
  location <- gd.list %>%
    filter(id =="1-s8aJtk6eUYqdAnqkYnYT-ryA50StDrq") %>%
    select(id)
  
  file_id <- as.character(location$id)
  out_path <- file.path(project_dir, "Arctic_PRISM_Tier1_plots.kmz")
  
  GET(url = paste0("https://www.googleapis.com/drive/v3/files/", file_id, "?alt=media"),
    config(token = gargle::token_fetch(scopes = "https://www.googleapis.com/auth/drive.readonly")),
    write_disk(out_path, overwrite = TRUE)
  )
}
  
species_lu <- read_csv(file.path(project_dir, "Arctic_PRISM_Species_Code.csv"))
survey_tbl <- read_csv(file.path(project_dir, "Arctic_PRISM_Surveys_Tiers1_Birds_Data_Donnees.csv"))

# PLOT (coodinates are in EPSG:3159. Need to be converted in EPSG:4326)
unzip(out_path, exdir = project_dir)
kml_file <- list.files(project_dir, pattern = "\\.kml$", full.names = TRUE)

xy <- st_read(kml_file) %>%dplyr::select(-Description)

xy <- xy %>%
  mutate(
    geometry = st_point_on_surface(geometry)
  ) %>%
  mutate(
    coords = st_coordinates(geometry)
  ) %>%
  mutate(
    lon = coords[,1],
    lat = coords[,2]
  ) %>%
  select(-coords) %>%
  st_drop_geometry()

 ## Validate mapping
 # library(ggplot2)
 # df_poly <- st_as_sf(xy, coords = c("lon", "lat"), crs = 4326)
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
  left_join(xy, by = c("Plot_Name" = "Name")) %>%
  filter(Species_Group != "Mammal", ECCC_species_code != "NOBI") %>%
  mutate(organization = "CWS",
       project= dataset_code,
       project_id	= NA,
       location	= paste0(dataset_code, ":", Plot_Name),
       location_id	= NA,
       location_buffer_m	= NA,
       latitude = lat,
       longitude	= lon,
       survey_id	= NA,
       survey_date	= as.Date(paste(Survey_year, sprintf("%02d", Survey_month), sprintf("%02d", Survey_day), sep = "-")),
       observer	= ifelse(is.na(ContactID), "obsNA", paste0("obs", abs(ContactID))),
       survey_distance_method	= "12 ha PRISM plot",
       survey_duration_method = "UNKNOWN",
       survey_duration	= if_else(!is.na(Finish_Time_1) & !is.na(Start_Time_1), as.character(as_hms(Finish_Time_1 - Start_Time_1)), NA_character_),
       detection_distance	= sub(" .*", "", Sighting_Type),
       detection_time	= "UNKNOWN",
       survey_time = as.character(format(strptime(Start_Time_1, "%H:%M"), "%H:%M:%S")),
       surveyDateTime = paste(survey_date, survey_time),
       survey_url	= NA,
       ind_count	= Total_Individuals_Adult + Total_Individuals_NonAdult,
       detection_heard	= "dnc",
       detection_seen	= "dnc",
       detection_comments = NA,
       #Behaviour
       originalBehaviourData = NA,
       pc_vt = "dnc",
       pc_vt_detail ="dnc",
       age = "dnc",
       fm = "dnc",
       group = "dnc",
       flyover = "dnc",
       displaytype = "dnc",
       nestevidence = "dnc",
       behaviourother = Count_Type
       )
  
  ################################
  #### Update master_observer ####
  ################################
  unique_observers <- survey %>%
    select(ContactID, observer) %>% 
    distinct() %>%
    filter(!is.na(ContactID)) %>% # Exclude rows where Observer is NA
    mutate(
      observer_name = as.character(ContactID),
      observer_id = observer
    )
  
  # Create the append_obs data frame
  append_obs <- unique_observers %>%
    select(observer_id, observer_name) %>%
    mutate(
      organization = "CWS-PAC",
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
  
  # Test species using scientific name
  species_tbl <- survey %>%
    select(ECCC_species_code) %>%
    dplyr::distinct(ECCC_species_code) %>%
    left_join(species_lu, by=c("ECCC_species_code"="2025_ECCC_Code")) %>%
    left_join(WT_spTbl, by=c("Scientific_Name"="scientific_name"))
  
  #list the ones that didn't pass
  species_tbl[is.na(species_tbl$species_code),]  
  
  survey <- survey %>% 
    left_join(species_tbl, by="ECCC_species_code")  %>%
    mutate(species = case_when(ECCC_species_code  == "AHGU" ~ "HERG", 
                               ECCC_species_code  == "CAGO_UNI" ~ "CANG",
                               ECCC_species_code  == "SHOR_UNI" ~ "UNSH",
                               ECCC_species_code  == "SCAU_UNI" ~ "UNSC", 
                               ECCC_species_code  == "PTAR_UNI" ~ "UNPT", 
                               ECCC_species_code  == "WATE_UNI" ~ "UNWT", 
                               ECCC_species_code  == "EIDE_UNI" ~ "UDIV", 
                               ECCC_species_code  == "GULL_UNI" ~ "UNGU", 
                               ECCC_species_code  == "JAEG_UNI" ~ "UJAE", 
                               ECCC_species_code  == "PASS_UNI" ~ "UNPA", 
                               ECCC_species_code  == "SPAR_UNI" ~ "UNSP", 
                               ECCC_species_code  == "BIRD_UNI" ~ "UNBI", 
                               ECCC_species_code  == "LOON_UNI" ~ "UNLO", 
                               ECCC_species_code  == "FALC_UNI" ~ "UNFA", 
                               ECCC_species_code  == "RAPT_UNI" ~ "URPT", 
                               ECCC_species_code  == "MEAD_UNI" ~ "UMED", 
                               ECCC_species_code  == "OWLS_UNI" ~ "UNOW", 
                               TRUE ~ species_code),
           species_scientific_name = Scientific_Name
           )

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
