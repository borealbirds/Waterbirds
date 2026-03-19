# Source dataset CWS Scoter Labrador Survey 2008-2009 (WLS2008-09)
# Author: "Melina Houle"
# Date: "February 10, 2026"
## Note on translation:
#     - Delete non bird obs (ARHA", "BEAV", "BLBE", "CARI", "MOOS", "MUSK", "PORC", "RFOX", "RIOT", "TIWO", "UNWH")
#     - 14 obs have count of 0 but birds was identified. Delete those. 
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
dataset_code <- "WLS2008-09"

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
dat <- "2008_09_enc_scoter_processed_d_ip.csv"
cond <- "2008_09_cond_scoter.csv"
sp <- "species_new.csv"

if (length(list.files(project_dir)) ==0) {
  gd.list <- drive_ls(as.character("https://drive.google.com/drive/folders/1HAeOJQK8K1ZKh54p9ZrJHxHKtQt7FSa_"))
  data <- gd.list %>%
    filter(name ==dat) %>%
    select("id")
  drive_download(as_id(as.character(data)), path = file.path(project_dir, dat))
  
  visit <- gd.list %>%
    filter(name ==cond) %>%
    select("id")
  drive_download(as_id(as.character(visit)), path = file.path(project_dir, cond))
  
  sp_lu <- gd.list %>%
    filter(name ==sp) %>%
    select("id")
  drive_download(as_id(as.character(sp_lu)), path = file.path(project_dir, sp))
}

sp_lu <- read.csv(file.path(project_dir, sp), fileEncoding="UTF-8-BOM")
visit <- read.csv(file.path(project_dir, cond), fileEncoding="UTF-8-BOM")
data <- read.csv(file.path(project_dir, dat), fileEncoding="UTF-8-BOM") %>%
  dplyr::filter(!is.na(X))

#xy_sf <- st_as_sf(data, coords = c("X", "Y"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()

# test on species
sp_valid <- data %>%
  left_join(sp_lu %>% select(spp_code, common_name,science_name),
            by = c("spp" = "spp_code")) %>%
  dplyr::filter(!spp %in% c("ARHA", "BEAV", "BLBE", "CARI", "MOOS", "MUSK", "PORC", "RFOX", "RIOT", "TIWO", "UNWH")) %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name),
            by = c("spp" = "species_code"))

invalid_codes <- sp_valid %>%
  distinct(spp, species_common_name) %>%
  filter(is.na(species_common_name)) 

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- data %>%
  dplyr::filter(!spp %in% c("ARHA", "BEAV", "BLBE", "CARI", "MOOS", "MUSK", "PORC", "RFOX", "RIOT", "TIWO", "UNWH")) %>%
  dplyr::filter(tot > 0) %>%
  left_join(visit %>% select(plot, yyyy, mm, dd, start1, end1, navigator), by = c("plot" = "plot", "year" = "yyyy")) %>%
  mutate(organization = "CWS-NOR",
         project= dataset_code,
         project_id	= NA,
         location	= paste0(organization, "_", dataset_code, ":", plot,":", index),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = Y,
         longitude	= X,
         survey_id	= NA,
         survey_time = ifelse(is.na(start1), "00:00:01", start1),
         survey_date	= paste0(year, "-", sprintf("%02d", as.integer(mm)), "-", sprintf("%02d", as.integer(dd)), " ", survey_time),
         survey_url	= NA,
         observer	= navigator,
         protocol_type = case_when(
           plot %% 1 == 0   ~ "AreaSearch_5x5",   
           plot %% 1 == 0.1 ~ "AreaSearch_2x2", 
           TRUE ~ NA_character_),
         survey_distance_method	= "UNKNOWN",
         survey_duration_method	= "UNKNOWN",
         start_time = strptime(start1, "%H:%M"),
         end_time = strptime(end1, "%H:%M"),
         dur_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
         survey_duration	= as.character(sprintf("%02d:%02d:%02d",
                                   dur_sec %/% 3600,
                                   (dur_sec %% 3600) %/% 60,
                                   dur_sec %% 60)),
         detection_distance	= "UNKNOWN",
         detection_time	= "UNKNOWN",
         species_code	= case_when(spp == "AGWT"  ~ "GWTE",
                                  spp == "UTER"  ~ "UNTE",
                                  spp == "CAGO"  ~ "CANG",
                                  spp == "TRSW"  ~ "UNBI",
                                  spp == "UNDI"  ~ "UDIV",
                                  spp == "USCA"  ~ "UNSC",
                                  spp == "HARD" ~ "HADU",
                                  spp == "MEGO" ~ "UNBI",
                                  spp == "NOGH" ~ "UNBI",
                                  spp == "UNPH" ~ "UPHL",
                                  TRUE ~ spp),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= tot,
         detection_heard	= "f",
         detection_seen	= "t",
         detection_comments = paste0(as.character(mal), " male- ", as.character(fem), " female- ", as.character(unk), " unknown")
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
WLS200809 <- survey %>%
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
save(WLS200809, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
