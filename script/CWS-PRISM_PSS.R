# Source dataset CWS PRISM WAPSUK National Park (ECCC)
# Author: "Melina Houle"
# Date: "February 2, 2025"
## Note on translation:
#         - survey follows PRISM protocol (plot survey).
#         - Walking Survey on 12 ha plot. 2 observer walking 25m apart on 400m transect. 
#         - Project occur in Manitoba (UTM 15N, EPSG:3159)
#         - Some obs seems duplicates but they aren't
#         - Instead of providing latitude/longitude, output stores lat long for eatch of the 4 corner of the plot. 
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
library(ggplot2)

# config.R store Google Drive ccredential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

#species
WT_spTbl <- read.csv(file.path("./lookupTables/species_codes.csv"),fileEncoding = "UTF-8-BOM") %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL"))

NC_spTbl <- read_excel(file.path("./lookupTables/lu_NatureCount_species.xlsx")) 

WT_durMethTbl <- read.csv(file.path("./lookupTables/duration_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distMethTbl <- read.csv(file.path("./lookupTables/distance_method_codes.csv"), fileEncoding="UTF-8-BOM")
WT_durBandTbl <- read.csv(file.path("./lookupTables/duration_interval_codes.csv"), fileEncoding="UTF-8-BOM")
WT_distBandTbl <- read.csv(file.path("./lookupTables/distance_band_codes.csv"), fileEncoding="UTF-8-BOM")

organization <- "CWS"
dataset_code <- "PRISM_PSS"
lu <- "./lookupTables"
project <- file.path("./project", dataset_code)

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
if (length(list.files(project)) ==0) {
  pid_pss <- "https://drive.google.com/drive/folders/1WO-tVYJM5fDyXqOX9OGPwVNd-ktce7Oq"
  gd.list <- drive_ls(as.character(pid_pss))
  survey_pss <- gd.list %>%
    filter(name =="prism-pss_naturecounts_data.txt") %>%
    select("id")
  drive_download(as_id(as.character(survey_pss)), path = file.path(project_dir, "pss_naturecounts_data.txt"))
}


pss <- read_delim(file.path(project_dir, "pss_naturecounts_data.txt"), delim = "\t", col_types = cols()) %>% 
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount))
  select(CollectionCode, ScientificName, DecimalLatitude, DecimalLongitude, UTMNorthing, UTMEasting, latitude, longitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
         Collector, SamplingEventIdentifier, SurveyAreaIdentifier, TimeObservationsStarted, EffortMeasurement2, ObservationCount, ObservationCount2, ObservationDescriptor2, ObservationCount3,
          ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, UTMZone,	UTMNorthing,	UTMEasting, CommonName, SpeciesCode)

#xy_sf <- st_as_sf(accws, coords = c("DecimalLongitude", "DecimalLatitude"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()


# test on species
accws_valid <- accws %>%
  distinct(SpeciesCode) %>%
  left_join(NC_spTbl %>% select(species_code, species_id),
            by = c("SpeciesCode" = "species_id"))
invalid_codes <- accws_valid %>%
  filter(is.na(species_code)) %>%
  distinct(SpeciesCode)

invalid_codes

unique(accws[accws$SpeciesCode %in% c(40432, 41519, 40419, 5950, 41200),
                    c("SpeciesCode", "CommonName")])
#################################################
## HARMONIZE
#################################################
survey <- accws %>%
  left_join(NC_spTbl %>% select(species_code, species_id),
            by = c("SpeciesCode" = "species_id")) %>%
  mutate(organization = "CWS",
         project= dataset_code,
         project_id	= NA,
         location	= paste0(dataset_code, ":", SurveyAreaIdentifier),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = DecimalLatitude,
         longitude	= DecimalLongitude,
         survey_id	= NA,
         survey_time = format(as.POSIXct(as.numeric((TimeCollected[!is.na(TimeCollected)])) * 3600, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S"),
         survey_date	= paste0(YearCollected, "-", MonthCollected, "-", DayCollected, " ", survey_time),
         survey_url	= NA,
         observer	= Collector,
         survey_distance_method	= "0m-100m-250m-500m-1000m",
         survey_duration	= NA,
         detection_distance	= case_when(EffortMeasurement2 =="<100m" ~ "0m-100m",
                                        EffortMeasurement2 =="100m" ~ "0m-100m",
                                        EffortMeasurement2 =="250m" ~ "0m-250m",
                                        EffortMeasurement2 =="500m" ~ "0m-500m",
                                        EffortMeasurement2 =="1 km" ~ "0m-1000m",
                                        TRUE ~ "UNKNOWN"),
         detection_time	= NA,
         species_code	= case_when(SpeciesCode == 40432  ~ "COMU",
                                  SpeciesCode == 41519  ~ "UNSH",
                                  SpeciesCode == 40419  ~ "UNTE",
                                  SpeciesCode == 5950  ~ "UNSH",
                                  SpeciesCode == 41200 ~ "UNME",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= ObservationCount,
         detection_heard	= "DNC",
         detection_seen	= "DNC",
         detection_comments = NA)# 







lf <- list.files(file.path(project_dir), full.names = TRUE)

all_data <- lf %>%
  lapply(read_delim, delim = "\t", col_types = cols(.default = col_character())) %>%  
  bind_rows() %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(is.na(NoObservations)) %>%
  select(CollectionCode, ScientificName, DecimalLatitude, DecimalLongitude, GeodeticDatum, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
       Collector, SamplingEventIdentifier, TimeObservationsStarted,TimeIntervalStarted, TimeIntervalEnded, ObservationCount, ObservationDescriptor,
       ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, ObservationCount5,  ObservationDescriptor5,
       ObservationCount6, ObservationDescriptor6, ObservationDate, UTMZone,	UTMNorthing,	UTMEasting, CommonName, Remarks2)



unique(all_data$ObservationDescriptor)







#---
WT.main.report <- survey %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude_NE,  longitude_NE, latitude_SE, longitude_SE,
                            latitude_SW, longitude_SW, latitude_NW , longitude_NW, survey_id, survey_date, survey_url, observer, survey_distance_method,
                            survey_duration, detection_distance, detection_time, survey_time, species_code, species_common_name, species_scientific_name,
                            ind_count, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "keep")

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
to_upload_contents <- drive_ls(as_id(dr)) # print(to_upload_contents)
upload_folder <- to_upload_contents[to_upload_contents$name == organization, ]
if (nrow(upload_folder) == 0) {
  upload_folder <- drive_mkdir(organization, path = as_id(dr))
}

# Save
save(WT.main.report, file = file.path(out_dir, paste0(dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(upload_folder), name = paste0(dataset_code,".RData"), overwrite = TRUE) 

# Create sub folder in 'toUpload' with the organization name
dr<- drive_get("Waterfowl/script", shared_drive = "BAM_AvianData")

script_path <- file.path(wd, "script", paste0(organization, "-", dataset_code,".R"))
drive_upload(media = script_path, path = as_id(dr), name = paste0(organization, "-", dataset_code,".R"), overwrite = TRUE) 
