# ---
# title: Waterbird models - Data - Wrangle Dataset ----
# author: Elly Knight
# created: March 8, 2026
# ---

#NOTES################################

#PURPOSE: This script combines the BAM dataset with waterbird-specific datasets for waterbird modelling.

#Some of the projects manually loaded as .Rdata objects will be uploaded to WildTrax and will be part of future versions of the BAM dataset, and so should be removed from this script

#This script currently only intakes point-survey data for the purpose of upland flyover removal and basic habitat relationships. Other data types will need to be added for future versions.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(wildrtrax) #WT data wrangling
library(avilistr) #taxonomy
library(data.table) #rbind wtih fill

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"
root_data <- "G:/Shared drives/BAM_AvianData"

#3. WT login ----
source("WTlogin.R")
wt_auth()

#GET DATA #####################

#1. BAM Dataaset ----
load(file.path(root_data, "BAMDataset", "04_BAMDataset_WT-2026-03-09_EBd-Jan-2026.Rdata"))

#2. BC MMP -----
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "BirdsCanada_BCMMP.RData"))
bcmmp <- WT.main.report
rm(WT.main.report)

#3. Yukon roadside survey ----
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "ECCC_YRBWS.RData"))
yrbws <- WT.main.report
rm(WT.main.report)

#4. Prairie pothole area search ----
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "USFWS_FSMS.RData"))
fsms <- WT.main.report
rm(WT.main.report)

#5. CWS aerial surveys ----
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "CWS-NOR_SES2017-19.RData"))
cwsnor1 <- WT.main.report
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "CWS-NOR_WLS1992.RData"))
cwsnor2 <- WT.main.report
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "CWS-NOR_WLS2008-09.RData"))
cwsnor3 <- WT.main.report
cwsnor <- rbind(cwsnor1, cwsnor2, cwsnor3)

#6. More aerial transects ----
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "USFWS_WBPHS.RData"))
wbphs <- WT.main.report
load(file.path(root_data, "Harmonization", "WildTrax-dataHarmonization", "waterfowl_Rdata", "ECCC_RLWS.RData"))
rlws <- WT.main.report
rm(WT.main.report)

#WATERBIRD SPECIES LIST ###############

#Using only species with 4-letter codes in WT

#1. Get WildTrax codes ----
#take out duplicates of scientific name
dup <- c("GRAJ", "CORBRA", "MEGU", "PICHUD", "ANSROS", "PSFL")

spp_wt <- wildrtrax::wt_get_species() |> 
  dplyr::filter(species_class=="AVES",
                species_scientific_name!=" ",
                !species_code %in% dup,
                str_length(species_code)==4,
                str_sub(species_code, 1, 2)!="45") |> 
  rename(scientific_name = species_scientific_name,
         common_name = species_common_name) |>  
  dplyr::select(species_code, common_name, scientific_name)

#2. Get avilist for families ----
fam <- c("Anatidae", "Scolopacidae", "Charadriidae", "Recurvirostridae", "Haematopodidae", "Phalarapodidae", "Rallidae", "Ardeidae", "Podicipedidae", "Gaviidae")
data(avilist_2025_short)
spp_avi <- avilist_2025_short |> 
  rename(scientific_name = Scientific_name) |> 
  dplyr::filter(Family %in% fam)

#3. Put together ----
spp <- spp_wt |> 
  inner_join(spp_avi) |> 
  inner_join(data.frame(species_code = colnames(dat))) |> 
  dplyr::select(all_of(c(colnames(spp_wt), "Family")))

#WRANGLE NEW DATASETS############

#1. BC MMP ----
# all species
bcmmp_wide <- bcmmp |> 
  wt_make_wide() |> 
  mutate(method = "PC - MMP",
         duration = as.integer(str_extract(survey_duration_method,
                                           "(?<=-)[0-9]+(?=min?)"))*60,
         distance = Inf,
         survey_date = ymd_hms(survey_date),
         project = "BCMMP",
         location_id = max(dat$location_id) + row_number()) |> 
  rename(date_time = survey_date) |> 
  dplyr::select(any_of(colnames(dat)))

#2. Yukon roadside survey ----
#waterfowl and 6 species as of 2004, we split to ensure those 6 get NAs
yrbws_wide_1 <- yrbws |> 
  dplyr::filter(year(ymd_hms(survey_date))< 2004) |> 
  mutate(survey_duration_method = NA) |> 
  wt_make_wide() |> 
  mutate(method = "PC",
         duration = NA,
         distance = Inf,
         survey_date = ymd_hms(survey_date),
         project = "YTRD_1",
         location_id = max(bcmmp_wide$location_id) + row_number()) |> 
  rename(date_time = survey_date) |> 
  dplyr::select(any_of(colnames(dat)))

yrbws_wide_2 <- yrbws |> 
  dplyr::filter(year(ymd_hms(survey_date)) >= 2004) |> 
  mutate(survey_duration_method = NA) |> 
  wt_make_wide() |> 
  mutate(method = "PC",
         duration = NA,
         distance = Inf,
         survey_date = ymd_hms(survey_date),
         project = "YTRD_2",
         location_id = max(yrbws_wide_1$location_id) + row_number()) |> 
  rename(date_time = survey_date) |> 
  dplyr::select(any_of(colnames(dat)))

#3. Aerial surveys ----
aerial_wide <- rbind(cwsnor, wbphs, rlws) |> 
  mutate(survey_duration_method = NA) |> 
  wt_make_wide() |> 
  mutate(method = "Aerial transect",
         duration = NA,
         distance = NA,
         survey_date = ymd_hms(survey_date),
         location_id = max(yrbws_wide_2$location_id) + row_number()) |> 
  rename(date_time = survey_date) |> 
  dplyr::select(any_of(colnames(dat)))

#4. Prairie pothole search ----
fsms_wide <- fsms |> 
  mutate(survey_duration_method = NA) |> 
  wt_make_wide() |> 
  mutate(method = "Plot search",
         duration = NA,
         distance = NA,
         survey_date = ymd_hms(survey_date),
         project = "BCMMP",
         location_id = max(aerial_wide$location_id) + row_number()) |> 
  rename(date_time = survey_date) |> 
  dplyr::select(any_of(colnames(dat)))
  
#5. Put together ----
all_1 <- rbindlist(list(dat, bcmmp_wide, yrbws_wide_2, aerial_wide, fsms_wide), fill=TRUE) |> 
  mutate(across(-colnames(dat[,1:10]), replace_na, 0))
all_2 <- rbindlist(list(all_1, yrbws_wide_1), fill=TRUE)

#6. Get just the waterbirbs ----
use <- all_2 |> 
  dplyr::select(any_of(c(colnames(dat[,1:10]), spp$species_code))) |> 
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude))

#7. Save -----
save(use, spp, file = file.path(root, "data", "01_WaterBirdData_Wide.Rdata"))

