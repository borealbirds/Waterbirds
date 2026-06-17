# Source dataset CWS Scoter Experimental Survey 1992 (WLS1992)
# Author: "Melina Houle"
# Date: "February 10, 2026"
## Note on translation:
#     - Delete red fox obs
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
dataset_code <- "WLS1992"

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
dat <- "WLS1992.csv"

if (length(list.files(project_dir)) ==0) {
  gd.list <- drive_ls(as.character("https://drive.google.com/drive/folders/1PSEgmuMYWRnDwsCaaAvWtD4J0_nXz1kg"))
  data <- gd.list %>%
    filter(name ==dat) %>%
    select("id")
  drive_download(as_id(as.character(data)), path = file.path(project_dir, dat))
}

data <- read.csv(file.path(project_dir, dat), fileEncoding="UTF-8-BOM") %>%
  dplyr::filter(lat_deg !=0)

#xy_sf <- st_as_sf(data, coords = c("lon_deg", "lat_deg"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()

# test on species
sp_valid <- data %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name),
            by = c("SPP" = "species_code"))
invalid_codes <- sp_valid %>%
  distinct(SPP, species_common_name) %>%
  filter(is.na(species_common_name)) 

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- data %>%
  dplyr::filter(SPP != "RFOX") %>%
  mutate(organization = "CWS-NOR",
         project= dataset_code,
         project_id	= NA,
         location	= paste0(organization, "_", dataset_code, ":", PlotOld,":", indx),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = lat_deg,
         longitude	= lon_deg,
         survey_id	= NA,
         survey_time = "00:00:01",
         survey_date	= paste0(YEAR, "-", sprintf("%02d", as.integer(MONTH)), "-", sprintf("%02d", as.integer(day)), " ", survey_time),
         survey_url	= NA,
         observer	= "Scott Gilliland",
         protocol_type = "AreaSearch_NA",
         survey_distance_method	= "UNKNOWN",
         survey_duration_method	= "UNKNOWN",
         survey_duration	= "UNKNOWN",
         detection_distance	= "UNKNOWN",
         detection_time	= "UNKNOWN",
         species_code	= case_when(SPP == "CAGO"  ~ "CANG",
                                  SPP == "AGWT"  ~ "GWTE",
                                  SPP == "USCA"  ~ "UNSC",
                                  SPP == "OLDS"  ~ "UNBI",
                                  SPP == "USCO"  ~ "USCT",
                                  SPP == "HARD" ~ "HADU",
                                  TRUE ~ SPP),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= MAL + FEM + UNK,
         detection_heard	= "f",
         detection_seen	= "t",
         detection_comments = paste0(as.character(MAL), " male- ", as.character(FEM), " female- ", as.character(UNK), " unknown")
  )
# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
WLS1992 <- survey %>%
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
save(WLS1992, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
