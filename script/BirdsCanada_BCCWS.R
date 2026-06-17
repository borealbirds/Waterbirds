# Source dataset Birds Canada BC Coastal Waterbird Surveys from Nature Count (BCCWS)
# Author: "Melina Houle"
# Date: "February 23, 2026"
## Note on translation:
#     - alcid sp. doesn't have WildTrax code. Translated as UDIV. species saved in comments
#     - Plover sp. doesn't have WildTrax code. Translated as UNSH. species saved in comments
#     - Eurasian x American Wigeon (hybrid). Translated as AMWI species saved in comments
#     - Mallard x Northern Pintail (hybrid). Translated as NOPI
#     - Data have a different structure than the ACCWS data also found on NatureCount (species and distance protocol)
#     - Species use acronym. Need to use search_species_code from the naturecounts R package
#    
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

# Wildtrax species
WT_spTbl <- read_csv(file.path("./lookupTables/species_codes.csv")) %>%
  dplyr::filter(!species_code %in% c("CORBRA", "PICHUD", "GRAJ", "PSFL", "STETRI", "183"))

# NatureCount species
#install.packages("naturecounts", 
#                 repos = c(birdscanada = 'https://birdscanada.r-universe.dev',
#                           CRAN = 'https://cloud.r-project.org'))
library(naturecounts)
sp_Tbl <- search_species_code(result= "all") 

organization <- "BirdsCanada"
dataset_code <- "BCCWS"

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
    pid_bccws <- "https://drive.google.com/drive/folders/1R5M2jj_aPSb59diGdBqkcAaAmkWROno0"
  gd.list <- drive_ls(as.character(pid_bccws))
  survey_bccws <- gd.list %>%
    filter(name =="bccws_naturecounts_data.txt") %>%
    select("id")
  drive_download(as_id(as.character(survey_bccws)), path = file.path(project_dir, "bccws_naturecounts_data.txt"))
}

bccws <- read_delim(file.path(project_dir, "bccws_naturecounts_data.txt"), delim = "\t", col_types = cols()) %>% 
  dplyr::select(where(~ !all(is.na(.)))) %>%
  dplyr::filter(!is.na(ObservationCount)) %>%
  select(CollectionCode, ScientificName, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected, DayCollected, TimeCollected,                             
       CollectorNumber, SurveyAreaIdentifier, TimeObservationsStarted, TimeObservationsEnded, EffortMeasurement2, ObservationCount, ObservationCount2, 
       ObservationDescriptor2, ObservationCount3,ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, CommonName, SpeciesCode)

# Validate XY
#xy_sf <- st_as_sf(bccws, coords = c("DecimalLongitude", "DecimalLatitude"))
#xy_sf <- st_set_crs(xy_sf, 4326)
#canada <- st_read("E:/MelinaStuff/BAM/GIS_layer/CanadaLAEA.shp")
#bnd <- st_transform(canada, crs= st_crs(4326))
#ggplot() +
#  geom_sf(data = bnd, fill = NA) +
#  geom_sf(data = xy_sf, fill = NA, colour = "red", linewidth = 6) +
#  theme_minimal()


# test on species
bccws_sp <- bccws %>%
  left_join(sp_Tbl %>% dplyr::select(BSCDATA, scientific_name, english_name), by = c("SpeciesCode" = "BSCDATA")) %>%
  mutate(english_name = case_when(CommonName == "Greater/Lesser Scaup" ~ "Unidentified Greater / Lesser Scaup",
                                  CommonName == "Iceland Gull (Thayer's)" ~ "Iceland Gull",
                                  CommonName == "Iceland Gull (glaucoides/kumlieni)" ~ "Iceland Gull",
                                  CommonName == "Pacific/Winter Wren" ~ "Unidentified Pacific / Winter Wren",
                                  CommonName == "Mallard (Domestic type)" ~ "Mallard",
                                  CommonName == "Larus sp." ~ "Unidentified Larus Gull",
                                  CommonName == "Common/Red-breasted Merganser" ~ "Unidentified Common/Red-breasted Merganser",
                                  CommonName == "peep sp." ~ "Unidentified Peeper Sandpiper",
                                  CommonName == "murrelet sp." ~ "Unidentified Murrelet" ,
                                  CommonName == "Cackling/Canada Goose" ~ "Cackling Goose",
                                  CommonName == "Merlin (Black)" ~ "Merlin",
                                  CommonName == "Herring x Glaucous-winged Gull (hybrid)" ~ "Western x Glaucous-winged Gull Hybrid",
                                  CommonName == "Herring Gull (American)" ~ "Herring Gull",
                                  CommonName == "Yellow-rumped Warbler (Audubon's)" ~ "Yellow-rumped Warbler",
                                  CommonName == "Trumpeter/Tundra Swan" ~ "Unidentified Swan",
                                  CommonName == "Tundra Swan (Whistling)" ~ "Tundra Swan",
                                  CommonName == "American Crow (American)" ~ "American Crow",
                                  CommonName == "Yellow-rumped Warbler (Myrtle)" ~ "Yellow-rumped Warbler",
                                  CommonName == "jaeger sp." ~ "Unidentified Jaeger",
                                  CommonName == "Common Snipe" ~ "Common Snipe",
                                  CommonName == "Peregrine Falcon" ~ "Peregrine Falcon",
                                  CommonName == "American/Northwestern Crow" ~ "American Crow",
                                  CommonName == "Dark-eyed Junco (Oregon)" ~ "Dark-eyed Junco",
                                  CommonName == "Short-billed Gull" ~ "Short-billed Gull",
                                  CommonName == "Mew Gull" ~ "Mew Gull",
                                  TRUE ~ english_name))

# test on species
bccws_valid <- bccws_sp %>%
  dplyr::distinct(SpeciesCode, scientific_name, english_name) %>%
  left_join(WT_spTbl, by = c("english_name"="species_common_name"))

invalid_codes <- bccws_valid %>%
  filter(is.na(species_code)) %>%
  distinct(english_name)

invalid_codes

#################################################
## HARMONIZE
#################################################
survey <- bccws_sp %>%
  left_join(WT_spTbl, by = c("english_name"="species_common_name")) %>%
  mutate(organization = "BirdsCanada",
         project= dataset_code,
         project_id	= NA,
         location = paste(paste0(organization,"_",dataset_code), SurveyAreaIdentifier, sep= ":"),
         location_id	= NA,
         location_buffer_m	= NA,
         latitude = DecimalLatitude,
         longitude = DecimalLongitude, 
         start_time = if_else(is.na(TimeObservationsStarted), as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"), as.POSIXct(as.numeric(TimeObservationsStarted) * 3600, origin = "1970-01-01", tz = "UTC")),
         end_time = if_else(is.na(TimeObservationsEnded), as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"), as.POSIXct(as.numeric(TimeObservationsEnded) * 3600, origin = "1970-01-01", tz = "UTC")),
         survey_id	= NA,
         survey_url	= NA,
         survey_time = format(start_time, "%H:%M:%S"),
         survey_date	= paste0(YearCollected, "-", sprintf("%02d", as.integer(MonthCollected)), "-", sprintf("%02d", as.integer(DayCollected)),  " ", survey_time),
         observer	= "obsNA",
         protocol_type = "Modified Transect",
         survey_distance_method	= "0m-100m-250m-500m-750m-1000m-INF",
         survey_duration_method	= "UNKNOWN",
         duration = round(as.numeric(difftime(end_time, start_time, units = "mins")),0),
         survey_duration	= case_when(duration <0 ~ "UNKNOWN",
                                     TRUE ~ paste0(as.character(duration),"min")),
         detection_distance	= case_when(EffortMeasurement2 == "<100m"~ "0m-100m",
                                        EffortMeasurement2 == "100m"~ "100m-250m",
                                        EffortMeasurement2 == "250m"~ "250m-500m",
                                        EffortMeasurement2 == "500m"~ "500m-750m",
                                        EffortMeasurement2 == "750m"~ "750m-1000m",
                                        EffortMeasurement2 == "1 km"~ "1000m-1250m",
                                        EffortMeasurement2 == "Unlimited"~ "1250m-INF",
                                        is.na(EffortMeasurement2) ~ "UNKNOWN"),
         detection_time	= "UNKNOWN",
         species_code	= case_when(english_name == "American Crow (Northwestern)"  ~ "AMCR",
                                  english_name == "gull sp."  ~ "UNGU",
                                  english_name == "loon sp." ~ "UNLO",
                                  english_name == "Shorebird sp."  ~ "UNSH",
                                  english_name == "scoter sp."  ~ "USCT",
                                  english_name == "tern sp."  ~ "UNTE",
                                  english_name == "Western x Glaucous-winged Gull (hybrid)"  ~ "GULL",
                                  english_name == "alcid sp."  ~ "UDIV",
                                  english_name == "duck sp."  ~ "UNDU",
                                  english_name == "Accipitrine hawk sp. (former Accipiter sp.)"  ~ "UAHA",
                                  english_name == "Phalarope sp."  ~ "UPHL",
                                  english_name == "merganser sp."  ~ "UNME",
                                  english_name == "Common/Barrow's Goldeneye" ~ "UGOL",
                                  english_name == "hawk sp." ~ "UAHA",
                                  english_name == "Bald Eagle"  ~ "BAEA",
                                  english_name == "American Herring Gull"  ~ "HERG",
                                  english_name == "cormorant sp."  ~ "UCOR",
                                  english_name == "Anas sp."  ~ "UNDU",
                                  english_name == "goose sp." ~ "UNGO",
                                  english_name == "Rock Pigeon (Feral Pigeon)" ~ "ROPI",
                                  english_name == "grebe sp."  ~ "UGRB",
                                  english_name == "Aythya sp."  ~ "UDIV",
                                  english_name == "Lesser/Greater Yellowlegs" ~ "UNYE",
                                  english_name == "Northern Pygmy-Owl"  ~ "NOPO",
                                  english_name == "Short-billed/Long-billed Dowitcher"  ~ "USLD",
                                  english_name == "bird sp."  ~ "UNBI",
                                  english_name == "Plover sp."  ~ "UNSH",
                                  english_name == "swan sp."  ~ "USWN",
                                  english_name == "falcon sp."  ~ "UNFA",
                                  english_name == "Eurasian x American Wigeon (hybrid)"  ~ "AMWI",
                                  english_name == "Hudsonian Whimbrel"  ~ "WHIM",
                                  english_name == "Green-winged Teal (Eurasian)"  ~ "GWTE",
                                  english_name == "Mallard x Northern Pintail (hybrid)" ~ "NOPI",
                                  english_name == "owl sp." ~ "UNOW",
                                  english_name == "Northern Yellow Warbler"  ~ "YEWA",
                                  english_name == "Graylag Goose (Domestic type)"  ~ "GRGO",
                                  english_name == "Buteo sp."  ~ "UNBU",
                                  english_name == "American/Pacific Golden-Plover (Lesser Golden-Plover)"  ~ "PAGP",
                                  english_name == "new world sparrow sp." ~ "UNSP",
                                  english_name == "Unidentified Murrelet" ~ "UNBI",
                                  TRUE ~ species_code),
         species_common_name	= WT_spTbl$species_common_name[match(species_code, WT_spTbl$species_code)],
         species_scientific_name	= WT_spTbl$scientific_name[match(species_code, WT_spTbl$species_code)],
         ind_count	= as.numeric(ObservationCount),
         detection_heard	= "dnc",
         detection_seen	= "t",
         detection_comments = ifelse(english_name %in% c("alcid sp.", "Plover sp.", "Eurasian x American Wigeon (hybrid)", "Mallard x Northern Pintail (hybrid)", "Unidentified Murrelet"), paste0("Original species: ", english_name), NA)
  )

# check
print(unique(survey$species_code[!(survey$species_code %in% WT_spTbl$species_code)]))

#---
BCCWS <- survey %>%
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
save(BCCWS, file = file.path(out_dir, paste0(organization,"_", dataset_code,".RData")))
report_path <- file.path(out_dir, paste0(organization,"_",dataset_code,".RData"))
drive_upload(media = report_path, path = as_id(dr), name = paste0(organization,"_",dataset_code,".RData"), overwrite = TRUE) 

script_path <- file.path(wd, "script", paste0(organization,"_",dataset_code,".R"))
file.copy(script_path, git_dir, overwrite = TRUE) 
