# Source dataset Birds Canada Coastal waterbird surveys from Nature Count (Coastal)
# Author: "Melina Houle"
# Date: "February 23, 2026"
## Note on translation:
#     - shearwater sp, alcid sp. don't have WildTrax code. Translated as UNBI and UDIV. species saved in comments
#     - Data have a different structure than the BCCWS data also found on NatureCount (species and distance protocol)
#     - Species use code. Need to use search_species from the naturecounts R package
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

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

# NatureCount species
#install.packages("naturecounts", 
#                 repos = c(birdscanada = 'https://birdscanada.r-universe.dev',
#                           CRAN = 'https://cloud.r-project.org'))
library(naturecounts)
sp_Tbl <- search_species(show = "all") %>%
  filter(class== "Aves")

organization <- "BirdsCanada"
dataset_code <- "ACCWS"

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

  pid_accws <- "https://drive.google.com/drive/folders/1poK7B5vou95gA9uag8oA9p0Tr1SqJG82"
  gd.list <- drive_ls(as.character(pid_accws))
  survey_accws <- gd.list %>%
    filter(name =="accws_naturecounts_data.txt") %>%
    select("id")
  drive_download(as_id(as.character(survey_accws)), path = file.path(project_dir, "accws_naturecounts_data.txt"))
  
}

accws <- read_delim(file.path(project_dir, "accws_naturecounts_data.txt"), delim = "\t", col_types = cols()) %>% 
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount)) %>%
  select(CollectionCode, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
       Collector, SurveyAreaIdentifier, TimeObservationsStarted, TimeObservationsEnded, EffortMeasurement2, ObservationCount, ObservationCount2, 
       ObservationDescriptor2, ObservationCount3,ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, CommonName, SpeciesCode)

# Validate XY
#xy_sf <- st_as_sf(accws, coords = c("DecimalLongitude", "DecimalLatitude"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()


accws_sp <- accws %>%
  left_join(sp_Tbl %>% dplyr::select(species_id, scientific_name, english_name), by = c("SpeciesCode" = "species_id"))  

# test on species
accws_valid <- accws_sp %>%
  dplyr::distinct(SpeciesCode, scientific_name, english_name) %>%
  left_join(WT_spTbl, by = "scientific_name")

invalid_codes <- accws_valid %>%
  filter(is.na(species_code)) %>%
  distinct(english_name)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- accws_sp %>%
  left_join(WT_spTbl %>% select(species_code, scientific_name), by = "scientific_name") %>%
  mutate(organization = "BirdsCanada",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), SurveyAreaIdentifier, sep= ":"),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = DecimalLatitude,
         longitude = DecimalLongitude, 
         start_time = as.POSIXct(as.numeric((TimeObservationsStarted[!is.na(TimeObservationsStarted)])) * 3600, origin = "1970-01-01", tz = "UTC"),
         end_time = as.POSIXct(as.numeric((TimeObservationsEnded[!is.na(TimeObservationsEnded)])) * 3600, origin = "1970-01-01", tz = "UTC"),
         survey_id	= NA,
         survey_time = format(start_time, "%H:%M:%S"),
         survey_date	= paste0(YearCollected, "-", sprintf("%02d", as.integer(MonthCollected)), "-", sprintf("%02d", as.integer(DayCollected)),  " ", survey_time),
         survey_url	= NA,
         observer	= Collector,
         protocol_type = "Modified Transect",
         survey_distance_method	= "0m-100m-250m-500m-1000m-INF",
         survey_duration_method	= "UNKNOWN",
         duration = round(as.numeric(difftime(end_time, start_time, units = "mins")),0),
         survey_duration	= case_when(duration <0 ~ "UNKNOWN",
                                     TRUE ~ paste0(as.character(duration),"min")),
         detection_distance	= case_when(EffortMeasurement2 == "<100m" ~ "0m-100m",
                                        EffortMeasurement2 == "100m"~ "100m-250m",
                                        EffortMeasurement2 == "250m"~ "250m-500m",
                                        EffortMeasurement2 == "500m"~ "500m-1000m",
                                        EffortMeasurement2 == "1 km"~ "1000m-INF",
                                        is.na(EffortMeasurement2) ~ "UNKNOWN"), 
         detection_time	= "UNKNOWN",
         species_code	= case_when(english_name == "Thick-billed/Common Murre"  ~ "UMUR",
                                  english_name == "shearwater sp."  ~ "UNBI",
                                  english_name == "gull sp." ~ "UNGU",
                                  english_name == "American Herring Gull"  ~ "HERG",
                                  english_name == "duck sp."  ~ "UNDU",
                                  english_name == "scoter sp."  ~ "USCT",
                                  english_name == "cormorant sp."  ~ "UCOR",
                                  english_name == "loon sp."  ~ "UNLO",
                                  english_name == "Common/Arctic Tern"  ~ "UNTE",
                                  english_name == "Shorebird sp."  ~ "UNSH",
                                  english_name == "alcid sp."  ~ "UDIV",
                                  english_name == "Greater/Lesser Scaup"  ~ "UNSC",
                                  english_name == "grebe sp." ~ "UGRB",
                                  english_name == "Common/Red-breasted Merganser" ~ "UNME",
                                  english_name == "tern sp."  ~ "UNTE",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= scientific_name,
         ind_count	= as.numeric(ObservationCount),
         detection_heard	= "f",
         detection_seen	= "t",
         detection_comments = ifelse(english_name %in% c("alcid sp.", "shearwater sp."), paste0("Original species: ", english_name), NA)
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
ACCWS <- survey %>%
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
save(ACCWS, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
