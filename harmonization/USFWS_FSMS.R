# Source dataset USWFS Four Square Mile Survey (FSMS)
# Author: "Melina Houle"
# Date: "February 16, 2026"
## Note on translation:
#     - The same pond can have multiple coordinates. They correspond to different pond within the plot. However, birds observed are link to pond number. 
#     - To allow the one to many relationship between pond coordinates and birds observed, I took the centroids of the pond network that shared the same unique id.
#     - More than 5500 ponds don't have coordinates, Use plot coordinates instead
#     - I used column protocol_type to report if location is Pond location or Plot location. 
#     - I reported Pond area in detection comments
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

organization <- "USFWS"
dataset_code <- "FSMS"

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
survey  <- "RawAllData.csv"
species <- "WaterfowlSpeciesLookupFSMS.csv"
pond <- "FsmsPondLookup.csv"
plot <- "FsmsPlotLookup.csv"

if (length(list.files(project_dir)) ==0) {
  
  survey_url <- "https://iris.fws.gov/APPS/ServCat/DownloadFile/275358"
  download.file(survey_url, file.path(project_dir, survey), mode = "wb")
  
  species_url <- "https://iris.fws.gov/APPS/ServCat/DownloadFile/270973"
  download.file(species_url, file.path(project_dir, species), mode = "wb")
  
  pond_url <- "https://iris.fws.gov/APPS/ServCat/DownloadFile/270975"
  download.file(pond_url, file.path(project_dir, pond), mode = "wb")
  
  plot_url <- "https://iris.fws.gov/APPS/ServCat/DownloadFile/270974"
  download.file(plot_url, file.path(project_dir, plot), mode = "wb")
}

survey_tbl <- read.csv(file.path(project_dir, survey), as.is = T, na.strings = c("NULL"), fileEncoding="UTF-8-BOM")
species_tbl <- read.csv(file.path(project_dir, species), colClasses = c(AOU = "character"), na.strings = c("NULL"), fileEncoding="UTF-8-BOM")
plot_tbl <- read.csv(file.path(project_dir, plot), as.is = T, na.strings = c("NULL"), fileEncoding="UTF-8-BOM") %>%
  dplyr::rename(XPlot = XCoord,
                YPlot = YCoord)

location_tbl <- read.csv(file.path(project_dir, pond), as.is = TRUE, na.strings = c("NULL", "NA"), fileEncoding = "UTF-8-BOM") %>%
  left_join(plot_tbl, by = c("PLOT" = "Plot")) %>%
  mutate(XCoord = as.numeric(XCoord),
         YCoord = as.numeric(YCoord),
         used_plot_coords = is.na(XCoord) | is.na(YCoord),
         x = coalesce(XCoord, XPlot), 
         y = coalesce(YCoord, YPlot)
  ) %>%
  group_by(Link) %>%
  summarise(
    # detect multiple coordinate pairs BEFORE averaging
    varied_coords = n_distinct(paste(x, y)) > 1,
    x = mean(x, na.rm = TRUE),
    y = mean(y, na.rm = TRUE),
    used_plot = any(used_plot_coords), # flag if any pond used plot coords
    pond_acres = mean(PondAcres, na.rm = TRUE), # summarize acres
    
    location_comments = case_when(
      used_plot ~ "Plot coordinates used because pond coordinates were missing",
      varied_coords ~ "Mean coordinates used due to multiple ponds",
      TRUE ~ paste0("Pond acres = ", round(pond_acres, 2))
    ),
    .groups = "drop"
  ) %>%
  # delete pond with no XY
  dplyr::filter(!is.na(x)) %>%
  # reproject in WGS 84
  st_as_sf(coords = c("x", "y"), crs = 5070) %>%
  st_transform(4326) %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude  = st_coordinates(.)[, 2]   
  ) %>%
  st_drop_geometry()

  ##Validate final XY
  #xy_sf <- st_as_sf(location_tbl, coords = c("longitude", "latitude"))
  #xy_sf <- st_set_crs(xy_sf, 4326)
  #usa <- st_read("E:/MelinaStuff/BAM/GIS_layer/NorthAmericaLazea.shp")
  #bnd <- st_transform(usa, crs= st_crs(4326))
  #ggplot() +
  #  geom_sf(data = bnd, fill = NA) +
  #  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
  #  theme_minimal()

fsms <- survey_tbl %>%
  inner_join(location_tbl, by = "Link") %>%
  left_join(species_tbl, by = "AOU" )


WT_spTbl_upper <- WT_spTbl %>%
  mutate(species_common_name = toupper(species_common_name))

# test on species
fsms_valid <- fsms %>%
  distinct(SPECIES) %>%
  left_join(WT_spTbl_upper %>% select(species_code, species_common_name),
            by = c("SPECIES" = "species_common_name"))
invalid_codes <- fsms_valid %>%
  filter(is.na(species_code)) %>%
  distinct(SPECIES)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- fsms %>%
  left_join(WT_spTbl_upper %>% select(species_code, species_common_name),
            by = c("SPECIES" = "species_common_name")) %>%
  mutate(organization = "USFWS",
         project= dataset_code,
         project_id	= NA,
         location	= paste0(dataset_code, ":", Link),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = latitude,
         longitude	= longitude,
         survey_id	= NA,
         survey_time = "00:00:01",
         survey_date	= format(as.Date(SurveyDate, format = "%m/%d/%Y"), "%Y-%m-%d"),
         survey_url	= NA,
         observer	= "obsNA",
         protocol_type = case_when(used_plot ~ "4 sq mi plot",
                                   TRUE ~ "Pond area survey"),
         survey_distance_method	= "UNKNOWN",
         survey_duration_method	= "UNKNOWN",
         survey_duration	= "UNKNOWN",
         detection_distance	= case_when(used_plot ~ "UNKNOWN",
                                        TRUE ~ "Pond"),
         detection_time	= "UNKNOWN",
         species_code	= case_when(SPECIES  == "NORHTERN PINTAIL"  ~ "NOPI",
                                  SPECIES  == "RED BREASTED MERGANSER" ~ "RBME",
                                  SPECIES  == "GREATER WHITE FRONTED GOOSE"  ~ "GWFG",
                                  is.na(SPECIES)  ~ "UNBI",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= Count,
         detection_heard	= "DNC",
         detection_seen	= "Yes",
         detection_comments = location_comments ) 

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
FSMS <- survey %>%
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
save(FSMS, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
