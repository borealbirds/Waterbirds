# Source dataset CWS Scoter Experimental Survey 2017-2019 (SES2017-19)
# Author: "Melina Houle"
# Date: "February 10, 2026"
## Note on translation:
#       - Delete all codes related to mammals: "OLDFUELCA", "AFOX", "BEARDEN", "BEAVER", "BLKBEAR", "CARIBOU", "GRIZ", "HARE", "MINK", "MOOSE", "MUOX", "OTTER,"OTTE", "PORCU", "UNKN", "WOLF", "WOLV"
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

organization <- "CWS-NOR"
dataset_code <- "SES2017-19"

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
dat <- "scoter_experimental_surveys_2017-2019_Data_Données.csv"
sp <- "scoter_experimental_surveys_2017-2019_SpeciesCodes_CodesDEspeces.csv"

if (length(list.files(project_dir)) ==0) {
  url_dat <- "https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fassess%2Fscoter-experimental-surveys_-_nt_-_mb_-_qc_-2017-2019%20%2Fscoter_experimental_surveys_2017-2019_Data_Donn%C3%A9es.csv"
  download.file(url_dat, file.path(project_dir, dat), mode = "wb")
  
  url_sp <- "https://data-donnees.az.ec.gc.ca/api/file?path=%2Fspecies%2Fassess%2Fscoter-experimental-surveys_-_nt_-_mb_-_qc_-2017-2019%20%2Fscoter_experimental_surveys_2017-2019_SpeciesCodes_CodesDEspeces.csv"
  download.file(url_sp, file.path(project_dir, sp), mode = "wb")
}

data <- read.csv(file.path(project_dir, dat), fileEncoding="UTF-8-BOM") %>%
  dplyr::filter(!is.na(Latitude))

sp_lu <- read.csv(file.path(project_dir, sp), fileEncoding="UTF-8-BOM")
#xy_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()

# test on species
sp_valid <- data %>%
  left_join(sp_lu %>% select(Project.code, X2024.Code, English.Name, Scientific.name), by = c("species" = "Project.code")) %>%
  dplyr::filter(!species %in% c("AFOX", "BEARDEN", "BEAVER", "BLKBEAR", "CARIBOU", "GRIZ", "HARE", "MINK", "MOOSE", "MUOX", "OTTER","OTTE", "PORCU", "UNKN", "WOLF", "WOLV")) %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name),
            by = c("X2024.Code" = "species_code"))
invalid_codes <- sp_valid %>%
  distinct(X2024.Code, species_common_name) %>%
  filter(is.na(species_common_name)) 

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- data %>%
  dplyr::filter(!species %in% c("OLDFUELCA","AFOX", "BEARDEN", "BEAVER", "BLKBEAR", "CARIBOU", "GRIZ", "HARE", "MINK", "MOOSE", "MUOX", "OTTER","OTTE", "PORCU", "UNKN", "WOLF", "WOLV")) %>%
  left_join(sp_lu %>% select(Project.code, X2024.Code, English.Name, Scientific.name), by = c("species" = "Project.code")) %>%
  mutate(organization = "CWS-NOR",
         project= dataset_code,
         project_id	= NA,
         location	= paste0(organization, "_", dataset_code, ":", studyarea,":", plot, ":", Feature_ID),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = Latitude,
         longitude	= Longitude,
         survey_id	= NA,
         survey_time = "00:00:01",
         survey_date	= paste0(y, "-", sprintf("%02d", as.integer(m)), "-", sprintf("%02d", as.integer(d)), " ", survey_time),
         survey_url	= NA,
         observer	= navigator,
         protocol_type = "AreaSearch_5x5",
         survey_distance_method	= "UNKNOWN",
         survey_duration_method	= "UNKNOWN",
         survey_duration	= "UNKNOWN",
         detection_distance	= "UNKNOWN",
         detection_time	= "UNKNOWN",
         species_code	= case_when(X2024.Code  == "CAGO"  ~ "CANG",
                                  X2024.Code  == "DIVE_UNI"  ~ "UDIV",
                                  X2024.Code  == "SCAU_UNI"  ~ "UNSC",
                                  X2024.Code  == "TERN_UNI"  ~ "UNTE",
                                  X2024.Code  == "DUCK_UNI"  ~ "UNDU",
                                  X2024.Code  == "GRLY_UNI"  ~ "UNYE",
                                  X2024.Code  == "PEEP_UNI" ~ "UPEE",
                                  X2024.Code  == "SCOT_UNI"  ~ "USCT",
                                  X2024.Code  == "MERG_UNI"  ~ "UNME",
                                  X2024.Code  == "SHOR_UNI"  ~ "UNSH",
                                  X2024.Code  == "UNGU_UNI" ~ "UNGU",
                                  X2024.Code  == "PTAR_UNI"  ~ "UNPT",
                                  X2024.Code  == "JAEG_UNI"  ~ "UJAE",
                                  X2024.Code  == "PHAL_UNI"  ~ "UPHL",
                                  X2024.Code  == "LOON_UNI"  ~ "UNLO",
                                  X2024.Code  == "GOLD_UNI"  ~ "UGOL",
                                  X2024.Code  == "DABB_UNI"  ~ "UDAB",
                                  X2024.Code  == "GREB_UNI"  ~ "UGRB",
                                  X2024.Code  == "SWAN"  ~ "USWN",
                                  X2024.Code  == "PLOV_UNI"  ~ "UNSH",
                                  X2024.Code  == "RAPT_UNI"  ~ "URPT",
                                  X2024.Code  == "REDB_UNI"  ~ "UNBI",
                                  is.na(X2024.Code) ~ "UNSH",
                                  TRUE ~ X2024.Code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= mal + fem + unk,
         detection_heard	= "f",
         detection_seen	= "t",
         detection_comments = paste0(as.character(mal), " male- ", as.character(fem), " female- ", as.character(unk), " unknown")
  )
# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
SES201719 <- survey %>%
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
save(SES201719, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
