# ---
# title: Waterbird models - Detectability - Habitat Summary
# author: Elly Knight
# created: March 13, 2026
# ---

#NOTES################################

#PURPOSE: This script uses two approaches to understand waterbird detection probability with ARUs and point counts ----

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(detect) #QPAD models
library(MuMIn) #AIC
library(data.table) #more data wrangling
library(dggridR) #spatial thinning
library(usdm)

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"

#3. Load waterbird data object ----
load(file.path(root, "data", "reporting_mar26", "04_WaterBirdData_NoFlyover.Rdata"))

options(scipen=9999)

#WRANGLING ###############

#1. Split into bird and covs, add temporal, standardize... ----
bird <- dat |> 
  dplyr::select(id, ABDU:YERA)

covs <- dat |> 
  dplyr::select(id, organization:distance) |> 
  mutate(hour = hour(date_time),
         day = yday(date_time)) |> 
  left_join(cov)

#2. Add the spatial grid ----
grid <- dgconstruct(spacing = 2.5, metric=TRUE)

covs$cell <- dgGEO_to_SEQNUM(grid, covs$longitude, covs$latitude)$seqnum

#3. VIF ----
covs_use <- covs |> 
  dplyr::select(uplandp_2km, occurrence_2km, recurrence_2km, seasonality_2km, wetp_2km)

vif(covs_use) #dang

#MODEL #############

#1. Set up a loop ----
spplist <- sort(colnames(bird |> dplyr::select(-id)))
for(i in 1:length(spplist)){
  
  spp.i <- spplist[i]
  
  #2. Get the data ----
  set.seed(i)
  cov.i <- covs |>
    group_by(cell) |>
    slice(sample.int(n(), 1)) |>
    ungroup()
  
  bird.i <- bird |> 
    dplyr::select(id, all_of(spp.i))
  colnames(bird.i) <- c("id", "count")
  
  dat.i <- inner_join(cov.i, bird.i) |> 
    mutate(occur = ifelse(count > 0, 1, 0)) |> 
    dplyr::filter(!is.na(occur),
                  !is.na(wetland),
                  !is.na(uplandp_2km),
                  !is.na(wetp_2km))
  
  #3. Check how many detections ----
  n <- sum(dat.i$occur)
  if(n < 200){next}
  
  #4. Model ----
  m.i <- glm(occur ~ wetland + uplandp_2km + occurrence_2km + seasonality_2km + wetp_2km + poly(hour, 2) + poly(day, 2), data = dat.i, family="binomial", na.action = "na.fail")
  
  #5. Save ----
  save(m.i, file=file.path(root, "results", "habitat_prelim", paste0(spp.i, ".Rdata")))
  
  cat(i, " ")
  
}

#SUMMARIZE############

#1. List of models ----
mods <- list.files(file.path(root, "results", "habitat_prelim"))

#2. Set up a loop ----
sum.list <- list()
for(i in 1:length(mods)){
  
  #3. Load model ----
  load(file.path(root, "results", "habitat_prelim", mods[i]))

  #4. Get summary ----
  sum.list[[i]] <- data.frame(summary(m.i)[["coefficients"]]) |> 
    rownames_to_column("covariate") |> 
    mutate(spp = str_sub(mods[i], 1, 4))
  
  cat(i, " ")
  
}

#5. Output ----
#rescale the betas zero to 1, remove nonsignificant ones
sum.out <- data.table::rbindlist(sum.list, fill=TRUE) |> 
  rename(p = 'Pr...z..', z = z.value, se = Std..Error) |> 
  mutate(Estimate = case_when(covariate=="uplandp_2km" ~ Estimate *100,
                              covariate=="wetp_2km" ~ Estimate*10,
                              !is.na(Estimate) ~ Estimate)) |> 
  mutate(across(c(Estimate, se, z, p), ~ round(.x, 3))) |> 
  rename(species_code = spp) |> 
  left_join(spp) |> 
  dplyr::filter(!covariate %in% c("(Intercept)", "poly(hour, 2)1", "poly(hour, 2)2", "poly(day, 2)1", "poly(day, 2)2"),
                p < 0.05) 

write.csv(sum.out, file.path(root, "results", "habitat_prelim_betas.csv"), row.names = FALSE)
