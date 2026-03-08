# ---
# title: Waterbird models - Detectability - Upland Truncation
# author: Elly Knight
# created: March 8, 2026
# ---

#NOTES################################

#PURPOSE: This script determines how to filter out upland surveys from the BAM dataset for waterbird modelling.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(terra) #rasters
library(sf) #shapefiles
library(wildrtrax) #species codes list
library(avilistr) #species wrangling
library(lme4) #mixed effects models
library(MuMIn) #dredging
library(flexmix) #Gaussian mixture models

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"
root_data <- "G:/Shared drives/BAM_AvianData"
root_covs <- "G:/Shared drives/BAM_NationalModels5/CovariateRasters"

#3. Load data object ----
load(file.path("BAMDataset", "04_BAMDataset_WT-2026-03-02_EBd-Jan-2026.Rdata"))

#4. Authenticate WildTrax ----
source("WTlogin.R")
wt_auth()

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
#Retain unidentifieds ----
unid <- c("UNCO", "UDOW", "UGOL", "UMAB", "UNGA", "UNLO", "UNSC", "UNWT", "UNYE", "UPHL", "USCT", "USLD", "UTEA", "UWCG")
spp <- spp_wt |> 
  inner_join(spp_avi) |> 
  dplyr::select(all_of(colnames(spp_wt))) |> 
  rbind(spp_wt |> 
          dplyr::filter(species_code %in% unid)) |> 
  inner_join(data.frame(species_code = colnames(dat))) 


#GET FLYOVER DATA ###############

#This section could probably be less jenky....

#1. GEt the WildTrax data ----
load(file.path(root_data, "BAMDataset", "WildTrax", "2026-03-02", "01_wildtrax_raw_2026-03-02.Rdata"))
pc <- do.call(rbind, pc.wt)

#2. Filter it ----
wt_fo <- pc |>  
  dplyr::filter(str_detect(detection_comments, "flyover"),
                species_code %in% spp$species_code)

#3. Get the IMBCR data ----
imbcr <- read.csv(file.path(root_data, "OtherDatasets", "IMBCR", "IMBCR_wrangled.csv")) |> 
  dplyr::filter(!How %in% c("-1", "-9", "U")) 

#4. Filter it ----
imbcr_fo <- imbcr |> 
  dplyr::filter(species %in% spp$species_code,
                How=="F") |> 
  rename(species_code = species)

#5. Put together -----
fo <- rbind(imbcr_fo |> 
              mutate(source = "IMBCR"),
            wt_fo |> 
              rename(datetime = survey_date,
                     lat = latitude,
                     lon = longitude,
                     distanceMethod = survey_distance_method,
                     durationMethod = survey_duration_method,
                     durationInterval = detection_time,
                     distanceBand = detection_distance,
                     count = individual_count) |> 
              mutate(How = "F",
                     Sex = "U",
                     id = row_number()) |> 
              dplyr::select(all_of(colnames(imbcr_fo))) |> 
              mutate(source = "WT"))

#6. Get the species that we should use ----
table(fo$species_code)

#7. Go get the non-flyover detections -----
wt_nofo <- pc |> 
  dplyr::filter(project_id %in% wt_fo$project_id,
                species_code %in% fo$species_code)

imbcr_nofo <- imbcr |> 
  dplyr::filter(species %in% fo$species_code,
                How!="F") |> 
  rename(species_code = species)

#8. Put together again ----
all <- rbind(imbcr_nofo |> 
               mutate(source = "IMBCR"),
             wt_nofo |> 
               rename(datetime = survey_date,
                      lat = latitude,
                      lon = longitude,
                      distanceMethod = survey_distance_method,
                      durationMethod = survey_duration_method,
                      durationInterval = detection_time,
                      distanceBand = detection_distance,
                      count = individual_count) |> 
               mutate(How = "O",
                      Sex = "U",
                      id = row_number()) |> 
               dplyr::select(all_of(colnames(imbcr_fo))) |> 
               mutate(source = "WT")) |> 
  mutate(flyover = 0) |> 
  rbind(fo |> 
          mutate(flyover = 1))

#9. Save out for GEE ----
write.csv(all, file.path(root, "data", "WaterbirdFlyovers.csv"), row.names = FALSE)

#MODEL ###################

#1. Get the covariates ----
cov_pt <- read.csv(file.path(root, "data", "WaterbirdFlyovers_water_point.csv")) |> 
  dplyr::select(-c(system.index, .geo, count)) |> 
  rename(dry_pt = b1,
         class_pt = b1_1,
         occur_pt = occurrence,
         season_pt = seasonality,
         recur_pt = recurrence) |> 
  dplyr::select(-c(b1_2, change_abs, change_norm, transition, max_extent))

cov_2km <- read.csv(file.path(root, "data", "WaterbirdFlyovers_water_20km.csv")) |> 
  dplyr::select(-c(system.index, .geo, count)) |> 
  rename(dry_2 = b1,
         class_2 = b1_1,
         occur_2 = occurrence,
         season_2 = seasonality,
         recur_2 = recurrence) |> 
  dplyr::select(-c(b1_2, change_abs, change_norm, transition, max_extent))

cov_20km <- read.csv(file.path(root, "data", "WaterbirdFlyovers_water_2km.csv")) |> 
  dplyr::select(-c(system.index, .geo, count)) |> 
  rename(dry_20 = b1,
         class_20 = b1_1,
         occur_20 = occurrence,
         season_20 = seasonality,
         recur_20 = recurrence) |> 
  dplyr::select(-c(b1_2, change_abs, change_norm, transition, max_extent))

#2. Put things together ----
all_cov <- all |> 
  inner_join(cov_pt) |> 
  inner_join(cov_200m) |> 
  inner_join(cov_2km) |> 
  mutate(across(c(occur_pt, occur_2, occur_20, recur_pt, recur_2, recur_20, season_pt, season_2, season_20), ~replace_na(.x, 0)),
         across(c(class_pt, class_2, class_20), ~as.factor(.x)))

#3. Data presents ----
all_wide <- all_cov |> 
  dplyr::select(flyover, dry_pt:season_20, -c(class_pt, class_2, class_20)) |> 
  pivot_longer(dry_pt:season_20, names_to="cov", values_to="val")

#gam
ggplot(all_wide) +
  geom_smooth(aes(x=val, y=flyover)) +
  facet_wrap(~cov, scales="free")

#linear
ggplot(all_wide) +
  geom_smooth(aes(x=val, y=flyover), method="lm") +
  facet_wrap(~cov, scales="free")

#Let's use:
#Linear: occur_2, season_2, occur_20, season_20, recur_20
#polynomial: dry_20

#4. Model no RE ----
glm1 <- glm(flyover ~ occur_2 + season_2 + occur_20 + season_20 + recur_20 + poly(dry_20),
              data=all_cov,
              family="binomial",
              na.action = "na.fail")

d1 <- dredge(glm1)
glm2 <- glm(flyover ~ occur_2 + season_2 + season_20 + recur_20,
            data=all_cov,
            family="binomial",
            na.action = "na.fail")
summary(glm2)

#5. Model with RE ----
glm3 <- glmer(flyover ~ occur_2 + season_2 + season_20 + recur_20 + (1|species_code),
              data=all_cov,
              family="binomial",
              na.action = "na.fail")
summary(glm3)

#6. Make predictions ----
all_pred <- data.frame(pred3 = predict(glm3, re.form = NA, type="response"),
                       pred3re = predict(glm3, type="response"),
                       pred2 = predict(glm2, type="response"),
                       pred = predict(glm1, type="response")) |> 
  cbind(all_cov)

wide_pred <- all_pred |> 
  dplyr::select(flyover, pred3, pred3re, pred2, pred) |> 
  pivot_longer(pred3:pred, names_to="model", values_to="prediction")

#7. Look at some things ----
ggplot(wide_pred) +
  geom_jitter(aes(x=prediction, y=flyover)) +
  facet_wrap(~model, ncol=2, scales="free")

#8. Use a Gaussian mixture model to find a truncation ----
set.seed(1234)
gmm <- flexmix(pred2 ~ 1 + species_code, data = all_pred, k=2, model = FLXMRglm(family = "gaussian"))
summary(gmm)

#9. Get the membership probabilities ----
all_pred$gmm_p <- posterior(gmm)[,1]
all_pred$gmm_fo <- clusters(gmm)

#10. Look at some things ----
table(all_pred$flyover, all_pred$gmm_fo)
