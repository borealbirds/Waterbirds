# ---
# title: Waterbird models - Detectability - Remove flyovers for upland areas
# author: Elly Knight
# created: March 14, 2026
# ---

#NOTES################################

#This script uses the output from 'detectablity/UplandFlyoverProbability.R' to remove likely flyovers from the overall waterbird dataset.

#PREAMBLE################

#1. Load packages----
library(tidyverse) #data wrangling
library(lme4) #mixed effects models
library(dismo) #BRTs
library(gbm) #BRT prediction
library(mgcv) #GAMs

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"

#3. Load truncation relations ----
load(file.path(root, "results", "UplandProbabilityOutput.Rdata"))

#4. Get the dataset ----
load(file.path(root, "data", "reporting_mar26", "02_WaterBirdData_Unknowns.Rdata"))
     
#WRANGLING ############

#1. Make long and join the spp info ----

#families that were modelled - let's predict to sandpipers for the oystercatchers
unique(all_pred$Family) 

dat <- clean |>
  pivot_longer(LBDO:CITE, names_to="species_code", values_to = "count") |>
  dplyr::filter(count > 0,
                !is.na(count)) |> 
  left_join(spp |> 
                dplyr::select(species_code, Family)) |> 
  mutate(Family = case_when(Family=="Haematopodidae" ~ "Scolopacidae",
                            !is.na(Family) ~ Family))

#2. Add the covariates ----
#rename for the BR
dat_cov <- dat |> 
  left_join(unique(cov)) |> 
  rename(class_pt = wc_class,
         dry_20 = uplandp_20km,
         occur_20 = occurrence_20km,
         season_20 = seasonality_20km,
         recur_20 = recurrence_20km,
         season_2 = seasonality_2km,
         dry_2 = uplandp_2km,
         occur_2 = occurrence_2km,
         recur_2 = recurrence_2km)
rm(dat)

#TRUNCATE THE DATASET ##############

#1. Make predictions ----
dat_cov$predbrt <- dismo::predict(object = brt, newdata=dat_cov, type="response")

#2. Set the threshold ----
ggplot(dat_cov) +
  geom_histogram(aes(x=predbrt)) +
  facet_wrap(~Family, scales="free")

thresh <- 0.2

#3. Try truncating ----
dat_cov$fo_pred <- ifelse(dat_cov$predbrt >= thresh, 1, 0)

#look at proportion
sum(dat_cov$fo_pred)/nrow(dat_cov)*100

#4. Remove from the dataset ----
removal <- dat_cov |> 
  mutate(count = ifelse(fo_pred==1, 0, count))

#5. Make it wide again ---- 
wide <- removal |> 
  dplyr::select(any_of(colnames(clean)), species_code, count) |> 
  arrange(species_code) |> 
  pivot_wider(names_from="species_code", values_from="count", values_fill=0) 


#6. Get the rest (that don't have waterbirds in them) ----
others <- dplyr::filter(clean, !id %in% wide$id)

#7. Join in back to the other data ----
dat <- data.table::rbindlist(list(wide, others), fill=TRUE) |> 
  arrange(id)

#7. Save ----
save(dat, cov, spp, file = file.path(root, "data", "reporting_mar26", "04_WaterBirdData_NoFlyover.Rdata"))
