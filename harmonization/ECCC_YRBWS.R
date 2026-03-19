# Source dataset ECCC Yukon Roadside Breeding Waterfowl Survey (1991-2016) (YRBWS)
# Author: "Melina Houle"
# Date: "February 18, 2026"
## Note on translation:
#     - 
#     - 
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

# config.R store Google Drive ccredential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

#species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

organization <- "ECCC"
dataset_code <- "YRBWS"

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
survey  <- "YT-RoadsideSurvey-1991-2016-Data-Donnée.csv"
species <- "YT-RoadsideSurvey-1991-2016-SpeciesList-ListDesEspeces.csv"

if (length(list.files(project_dir)) ==0) {

  survey_url <- "https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fassess%2Fyukon-co-operative-roadside-breeding-waterfowl-waterbird-survey-1991-2016%2FYT-RoadsideSurvey-1991-2016-Data-Donn%C3%A9e.csv"
  download.file(survey_url, file.path(project_dir, survey), mode = "wb")

  species_url <- "https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fassess%2Fyukon-co-operative-roadside-breeding-waterfowl-waterbird-survey-1991-2016%2FYT-RoadsideSurvey-1991-2016-SpeciesList-ListDesEspeces.csv"
  download.file(species_url, file.path(project_dir, species), mode = "wb")

}

survey_tbl <- read_csv(file.path(project_dir, survey)) %>%
  dplyr::filter(!spec_abbrev %in% c("NOOO", "ONZI", "OVDA"))

species_tbl <- read_csv(file.path(project_dir, species)) 

# Validate XY
#xy_sf <- survey_tbl %>%
#  distinct(loc_id, latitude, longitude) %>%
#  st_as_sf(coords = c("longitude", "latitude")) %>% st_set_crs(xy_s, 4326)
#bnd <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp") %>% st_transform(crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()

yrbws <- survey_tbl %>%
  left_join(species_tbl %>% dplyr::select(ProjectCode, EnglishName), by = c("spec_abbrev" = "ProjectCode")) %>%
  mutate(join_name = str_to_upper(str_trim(EnglishName)))

# test on species
yrbws_valid <- yrbws %>%
  distinct(spec_abbrev, EnglishName, join_name)%>%
  left_join(WT_spTbl %>% select(species_code, species_common_name) %>%
      mutate(join_name = str_to_upper(str_trim(species_common_name))), by = "join_name")

invalid_codes <- yrbws_valid %>%
  filter(is.na(species_code)) %>%
  distinct(spec_abbrev)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- yrbws %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name) %>%
              mutate(join_name = str_to_upper(str_trim(species_common_name))), by = "join_name") %>%
  mutate(organization = "ECCC",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), route_id, loc_id, sep= ":"),
         location_id	= NA,
         location_buffer_m	= NA,
         survey_id	= NA,
         survey_time = ifelse(hour_start==0, "00:00:01", as.character(sprintf("%02d:%02d:00", hour_start, minute_start))),
         survey_date	= paste0(year, "-", sprintf("%02d", as.integer(month)), "-", sprintf("%02d", as.integer(day)), " ", survey_time),
         survey_url	= NA,
         observer	= paste0("obs_", as.character(primary_obs)),
         protocol_type = "Wetland survey along the highway",
         survey_distance_method	= "UNKNOWN",
         survey_duration_method	= "UNKNOWN",
         duration = (hour_end * 60 + minute_end) - (hour_start * 60 + minute_start),
         survey_duration	= case_when(duration ==0 ~ "UNKNOWN",
                                     TRUE ~ as.character(duration)),
         detection_distance	= "UNKNOWN",
         detection_time	= "UNKNOWN",
         species_code	= case_when(spec_abbrev == "SCAU"  ~ "UNSC",
                                  spec_abbrev == "GOLD"  ~ "UGOL",
                                  spec_abbrev == "DOWI" ~ "UDOW",
                                  spec_abbrev == "MERG"  ~ "UNME",
                                  spec_abbrev == "SWAN"  ~ "USWN",
                                  spec_abbrev == "CAGO"  ~ "CANG",
                                  spec_abbrev == "YELL"  ~ "UNYE",
                                  spec_abbrev == "HBCI"  ~ "CITE",
                                  spec_abbrev == "BUTE"  ~ "UNHA",
                                  spec_abbrev == "CORE"  ~ "REDP",
                                  spec_abbrev == "PEEP"  ~ "UPEE",
                                  spec_abbrev == "PLOV"  ~ "UNSH",
                                  spec_abbrev == "EAGL" ~ "BAEA",
                                  spec_abbrev == "CAGS" ~ "CANG",
                                  spec_abbrev == "HYWI"  ~ "AMWI",
                                  spec_abbrev == "REDP"  ~ "UNRE",
                                  spec_abbrev == "TRIN"  ~ "UNYE",
                                  spec_abbrev == "WIGE"  ~ "AMWI",
                                  spec_abbrev == "GODW" ~ "UNSH",
                                  spec_abbrev == "CAGL" ~ "CANG",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= as.numeric(count),
         detection_heard	= "dnc",
         detection_seen	= "t",
         detection_comments = ifelse(spec_abbrev %in% c("WIGE", "HYWI", "PLOV", "GODW"), paste0("Original species: ", EnglishName), NA)
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
YRBWS <- survey %>%
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

# Save .Rdata
save(YRBWS, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

## Save WT upload location 
location <- YRBWS %>%
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
survey <- YRBWS %>%
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