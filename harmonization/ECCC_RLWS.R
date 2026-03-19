# Source dataset ECCC Rasmussen Lowlands Waterfowl survey, Nunavut 1994-1995 (RLWS)
# Author: "Melina Houle"
# Date: "February 18, 2026"
## Note on translation:
#     - CSV provided on ECCC datamart has a broken link, Need to use the zipped gdb
#     - Delete none bird obs ("AFOX", "none", "S_HARE", "WOLF")
#     - Protocol 2km transect 
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
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL"))

organization <- "ECCC"
dataset_code <- "RLWS"

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
survey  <- "Rasmussen_Lowlands_Waterfowl_Surveys_1994_1995_v2.gdb"
species <- "Rasmussen-Lowlands-Waterfowl-Surveys-1994-1995-SpeciesCodes-CodesDEspeces.csv"

if (length(list.files(project_dir)) ==0) {
  
  survey_url <- "https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fassess%2Fdistribution-and-abundance-of-waterfowl-in-the-rasmussen-lowlands-nunavut-1994-1995%2FRasmussen_Lowlands_Waterfowl_Surveys_1994_1995_v2.gdb.zip"
  download.file(survey_url, file.path(project_dir, paste0(survey, ".zip")), mode = "wb")
  unzip(file.path(project_dir, paste0(survey, ".zip")), exdir = project_dir)  

  species_url <- "https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fassess%2Fdistribution-and-abundance-of-waterfowl-in-the-rasmussen-lowlands-nunavut-1994-1995%2FRasmussen-Lowlands-Waterfowl-Surveys-1994-1995-SpeciesCodes-CodesDEspeces.csv"
  download.file(species_url, file.path(project_dir, species), mode = "wb")

}

survey_sf <- st_read(file.path(project_dir, survey), layer = "Rasmussen_Lowlands_Waterfowl_Surveys_1994_1995_v2") %>%
  dplyr::filter(!SPECIES %in% c("AFOX", "None", "S_HARE", "WOLF"))

species_tbl <- read_csv(file.path(project_dir, species))

  ##Validate  XY
  #xy_sf <- survey_sf 
  #usa <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp")
  # bnd <- st_transform(usa, crs= st_crs(4326))
  #ggplot() +
  #  geom_sf(data = bnd, fill = NA) +
  #  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
  #  theme_minimal()

rlws <- survey_sf %>%
  st_drop_geometry() %>%
  left_join(species_tbl %>% dplyr::select(ProjectCode, EnglishName), by = c("SPECIES" = "ProjectCode"))

# test on species
rlws_valid <- rlws %>%
  distinct(SPECIES, EnglishName) %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name),
            by = c("EnglishName" = "species_common_name"))
invalid_codes <- rlws_valid %>%
  filter(is.na(species_code)) %>%
  distinct(SPECIES)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- rlws %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name),
            by = c("EnglishName" = "species_common_name")) %>%
  mutate(organization = "ECCC",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), TRANSECT, SEGMENT, sep= ":"),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = LAT_DD_WGS84,
         longitude	= LONG_DD_WGS84,
         survey_id	= NA,
         survey_time = "00:00:01",
         survey_date	= paste0(SURVEY_YR, "-", sprintf("%02d", as.integer(SURVEY_MONTH)), "-", sprintf("%02d", as.integer(SURVEY_DAY)), " ", survey_time),
         survey_url	= NA,
         observer	= NA,
         protocol_type = "Area survey 2km transect",
         survey_distance_method	= "0m-1000m",
         survey_duration_method	= "UNKNOWN",
         survey_duration	= "UNKNOWN",
         detection_distance	= "UNKNOWN",
         detection_time	= "UNKNOWN",
         species_code	= case_when(SPECIES == "GEESE"  ~ "UNGO",
                                  SPECIES == "GEES"  ~ "UNGO",
                                  SPECIES == "CAN" ~ "CANG",
                                  SPECIES == "SNGOWH"  ~ "LSGW",
                                  SPECIES == "SNGOBP"  ~ "SNGO",
                                  SPECIES == "BIRD"  ~ "UNBI",
                                  SPECIES == "RAPTOR"  ~ "URPT",
                                  SPECIES == "SNGO"  ~ "SNGO",
                                  SPECIES == "GULL"  ~ "UNGU",
                                  SPECIES == "EIDER"  ~ "UDIV",
                                  SPECIES == "DIVER"  ~ "UDIV",
                                  SPECIES == "DUCKS"  ~ "UNDU",
                                  SPECIES == "SEADUCK" ~ "UDIV",
                                  SPECIES == "PTAR" ~ "UNPT",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= as.numeric(COUNT),
         detection_heard	= "dnc",
         detection_seen	= "t",
         detection_comments = case_when(SPECIES == "SNGOBP" ~ "blue morph",
                                        SPECIES == "EIDER" ~ "unidentified eider",
                                        SPECIES == "SEADUCK" ~ "unidentified seaduck")
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
RLWS <- survey %>%
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
save(RLWS, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
