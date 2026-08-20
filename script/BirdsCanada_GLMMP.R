# Source dataset Birds Canada Great Lakes Marsh Monitoring Program (GLMMP)
# Author: "Melina Houle"
# Date: "August 20, 2026"
## Note on translation:
#    -glmmp has 2 types of protocol: 0-5-10min and 0-5-10-15min. project would need to be split if stored in Wildtrax. .
#    -ObservationCount ="Total Count"
#    -ObservationCount2 ="0-5 min"
#    -ObservationCount3 ="5-10 min"
#    -ObservationCount4 ="10-15 min"
#    -ObservationCount5 ="Target species (5-15 min)"
#    -ObservationCount6 ="Outside/Flythrough"
#    - Some obs have total count withou any band. Need to be flag before pivot
#    - drop 18000 non-birds obs
#    - Distance band: 100m according to methodology
#    - Duration ban 0-5min used broadcast
#    - Duration band 5-10min and 10-15min were silent. 
#  Based on the description of ObservationCount columns, we used totalCount for abundance and try to split by time interval 0-5, 5-10 or 10-15 when known. 
#  We log the flyover but 
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
library(tibblify)
library(reshape2)

install.packages("naturecounts", 
                 repos = c(birdscanada = 'https://birdscanada.r-universe.dev',
                           CRAN = 'https://cloud.r-project.org'))
library(naturecounts)

# config.R store Google Drive credential and path to working directory
source("./config.R")

## Initialize variables (wd is define in config.R)
setwd(file.path(wd))

drive_auth()

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

# NatudeCount species
sp_Tbl <- read_excel(file.path("./lookupTables/lu_NatureCounts_species.xlsx"))

# Load NatureCount species list and drop all non birds
nc_sp_code <- search_species() %>% dplyr::filter(taxon_group== "BIRDS")

organization <- "BirdsCanada"
dataset_code <- "GLMMP"

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
  glmmp <- nc_data_dl(collections = "MMPBIRDS", fields_set = "extended", request_id = r_id, username = natureCounts_username, info = "download Atlas")
  save(glmmp, file = file.path(project_dir, "mmpbirds.RData"))
  #pid_glmmp <- "https://drive.google.com/drive/u/0/folders/1W0hCa4kgnuNobNa-xFfkc4afP9ABD22d"
  #gd.list <- drive_ls(as.character(pid_glmmp))
  #survey_glmmp <- gd.list %>%
  #  filter(name =="GLMMP.Rdata") %>%
  #  select("id")
  #drive_download(as_id(as.character(survey_glmmp)), path = file.path(project_dir, "GLMMP.Rdata"))
}else{
  load(file.path(project_dir, "mmpbirds.RData"))
}

#load(file.path(project_dir, "mmpbirds.Rdata"))
#write.csv(glmmp, "./out/BirdsCanada_GLMMP/raw.csv")

survey <- glmmp %>% 
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount)) %>%
  dplyr::filter(!is.na(DecimalLatitude)) %>%
  dplyr::inner_join(sp_Tbl, by = "species_id")  %>%# drop all non-birds
  dplyr::select(GlobalUniqueIdentifier, SamplingEventIdentifier, CollectionCode, latitude, longitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
                CollectorNumber, SurveyAreaIdentifier, TimeObservationsStarted, EffortMeasurement1, ObservationCount, ObservationCount2, 
                ObservationCount3, ObservationCount4, ObservationCount5, ObservationCount6, ScientificName, CommonName, SpeciesCode)

# Validate XY
library(ggplot2)
xy_sf <- st_as_sf(survey, coords = c("longitude", "latitude"))
xy_sf <- st_set_crs(xy_sf, 4326)
canada <- st_read("./lookupTables/GIS_layers/CanadaLAEA.shp")
bnd <- st_transform(canada, crs= st_crs(4326))
ggplot() +
  geom_sf(data = bnd, fill = NA) +
  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
theme_minimal()

# test on species
glmmp_valid <- survey %>%
  dplyr::distinct(SpeciesCode, CommonName) %>%
  rename(species_code = SpeciesCode) %>%
  left_join(WT_spTbl, by = "species_code")

invalid_codes <- glmmp_valid %>%
  filter(is.na(CommonName)) %>%
  distinct(species_code)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- survey %>%
  mutate(ObservationDescription7 = "unknown band",
         ObservationCount7 = as.character(as.numeric(ObservationCount) -rowSums(across(c(ObservationCount2, ObservationCount3,
                                                                                         ObservationCount4),
                                                                                       as.numeric), na.rm = TRUE)
                                          ))
  
survey_expanded <- melt(survey, measure.vars = c("ObservationCount2","ObservationCount3","ObservationCount4","ObservationCount6", "ObservationCount7"), value.name = "abundance")
survey_expanded$abundance <- as.numeric(survey_expanded$abundance)
survey_expanded <- subset(survey_expanded, survey_expanded$abundance>0 )

survey.flat <- survey_expanded %>%
  rename(species_code = SpeciesCode) %>%
  left_join(WT_spTbl %>% select(species_code, species_common_name, scientific_name), by = "species_code") %>%
  mutate(organization = "BirdsCanada",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), SurveyAreaIdentifier, sep= ":"),
         location_id	= NA,
         location_buffer_m = NA,
         survey_id	= NA,
         survey_time = ifelse(is.na(TimeCollected), "00:00:01", format(as.POSIXct(as.numeric((TimeCollected[!is.na(TimeCollected)])) * 3600, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S")),
         survey_date	= paste0(YearCollected, "-", sprintf("%02d", as.integer(MonthCollected)), "-", sprintf("%02d", as.integer(DayCollected)),  " ", survey_time),
         survey_url	= NA,
         observer	= "obsNA",
         protocol_type = "Broadcast Point count",
         survey_distance_method	= "0m-100m",
         survey_duration	= case_when(EffortMeasurement1 == "10" ~ "10min",
                                            EffortMeasurement1 == "15" ~ "15min"),
         survey_duration_method	= case_when(EffortMeasurement1 == "10" ~ "0-5-10min",
                                            EffortMeasurement1 == "15" ~ "0-5-10-15min"),
         detection_distance	= "0m-100m", 
         detection_time	= case_when(variable == "ObservationCount2" ~ "0-5min",
                                    variable == "ObservationCount3" ~ "5-10min",
                                    variable == "ObservationCount4" ~ "10-15min",
                                    variable == "ObservationCount6" ~ "flyover",
                                    variable == "ObservationCount7" ~ "UNKNOWN"),
         species_code	= case_when(CommonName == "Northern Yellow Warbler" ~ "YEWA",
                                  CommonName == "Canada Goose" ~ "CANG",
                                  CommonName == "moorhen/coot/gallinule sp."  ~ "UNGA",
                                  CommonName == "American Herring Gull"  ~ "AHGU",
                                  CommonName == "American Bittern"  ~ "AMBI",
                                  CommonName == "Green-winged Teal (American)"  ~ "ANACRE",
                                  CommonName == "Leach's Storm-Petrel"  ~ "LESP" ,
                                  CommonName == "Ring-necked Pheasant"  ~ "RNEP" ,
                                  CommonName == "Rock Pigeon (Feral Pigeon)"  ~ "ROPI",
                                  CommonName == "Nelson's Sparrow"  ~  "NESP",
                                  CommonName == "Gray Partridge"  ~  "GRAP",
                                  CommonName == "Eastern Whip-poor-will" ~  "EWPW",
                                  CommonName == "Dark-eyed Junco" ~  "DEJU",
                                  CommonName == "Canada Jay" ~  "CAJA",
                                  CommonName == "Barred Owl" ~  "BAOW",
                                  CommonName == "swallow sp."  ~ "UNSW",
                                  CommonName == "Marsh Wren" ~  "MAWR",
                                  CommonName == "Redpoll (flammea)" ~ "UNRE",
                                  TRUE ~ species_code),
         species_scientific_name	= scientific_name,
         ind_count	= abundance,
         detection_heard	= "dnc",
         detection_seen	= "dnc",
         detection_comments = "Broadcast 0-5min, silent 5-15min"
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
#---
GLMMP <- survey.flat %>%
  dplyr::filter(!detection_time == "flyover") %>%
  dplyr::group_by(organization, project,  project_id, location, location_id,  location_buffer_m, latitude,  longitude, survey_id, survey_date, survey_url, 
                  observer, survey_distance_method, survey_duration_method, survey_duration, detection_distance, detection_time, survey_time, species_code, 
                  species_common_name, species_scientific_name, detection_heard, detection_seen, detection_comments) %>%
  dplyr::summarise(individual_count = sum(ind_count), .groups= "drop")

#--------------------------------------------------------------
#
#       EXPORT
#
#--------------------------------------------------------------

# Create sub folder in 'toUpload' with the organization name
dr<- drive_get("waterfowl_Rdata/", shared_drive = "BAM_AvianData")

# Save .Rdata
save(GLMMP, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 
