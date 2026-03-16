# ---
# title: Waterbird models - Detectability - Detectability
# author: Elly Knight
# created: March 13, 2026
# ---

#NOTES################################

#PURPOSE: This script uses two approaches to understand waterbird detection probability with ARUs and point counts ----

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(detect) #QPAD models

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"

#3. Load waterbird data object ----
load(file.path(root, "data", "reporting_mar26", "04_WaterBirdData_NoFlyover.Rdata"))

options(scipen=9999)

#WRANGLING ###############

#1. Split into bird and covs, add temporal ----
bird <- dat |> 
  dplyr::select(id, ABDU:YERA)

covs <- dat |> 
  dplyr::select(id, organization:distance) |> 
  mutate(sensor = ifelse(method %in% c("eBird", "PC", "PC - MMP"), "PC", "ARU"),
         year = year(date_time)) |> 
  left_join(cov)

#SUMMARIZE DATASET#############

#1. Set up to do list ----
loop <- expand.grid(spp = as.character(colnames(bird |> dplyr::select(-id))),
                    sensor = c("PC", "ARU"))
out <- data.frame()
for(i in 1:nrow(loop)){
  
  spp.i <- loop$spp[i]
  sensor.i <- loop$sensor[i]
  
  #2. Filter data ----
  dat.i <- covs |> 
    left_join(bird |> 
                dplyr::select(id, all_of(spp.i)) |> 
                setNames(c("id", "species"))) |> 
    dplyr::filter(sensor==sensor.i,
                  !is.na(species)) |> 
    mutate(occur = ifelse(species > 0, 1, 0))
  
  #3. Survey occurrence ----
  survey.occur <- sum(dat.i$occur)/nrow(dat.i)
  
  #3. Site occurrence ----
  site.i <- dat.i |> 
    group_by(location_id, year) |> 
    summarize(count = sum(species),
              surveys = n()) |> 
    ungroup() |> 
    mutate(occur = ifelse(count > 0, 1, 0))
  
  site.occur <- sum(site.i$occur)/nrow(site.i)
  
  #4. Mean detection rate at sites with detections ----
  detrate.i <- site.i |> 
    dplyr::filter(occur==1) |> 
    left_join(dat.i |> 
                dplyr::select(-occur)) |>  
    mutate(occur = ifelse(species > 0, 1, 0)) |> 
    group_by(location_id, year) |> 
    summarize(dets = sum(occur),
              surveys = n()) |> 
    ungroup() |> 
    mutate(det = dets/surveys) |> 
    summarize(rate_mn = mean(det),
              rate_sd = sd(det))
  
  #5. Mean count for detections ----
  det.i <- dat.i |> 
    dplyr::filter(species > 0)
  count.i <- mean(det.i$species)
  count.sd.i <- sd(dat.i$species)

  #6. Put together ----
  out <- rbind(out,
               data.frame(spp = spp.i,
                          sensor = sensor.i,
                          survey_occur = survey.occur,
                          site_occur = site.occur,
                          occur_count = count.i,
                          occur_count_sd = count.sd.i,
                          occur_det = detrate.i$rate_mn,
                          occur_det_sd = detrate.i$rate_sd))
  
  cat(i, " ")
  
  #7. Save ----
  write.csv(out, file.path(root, "results", "OccurrenceRates.csv"), row.names=FALSE)
  
}

#DATA PRESENTS ########

#1. Wrangle ----
wide <- out |> 
  dplyr::select(-occur_count_sd, -occur_det_sd) |> 
  pivot_longer(survey_occur:occur_det, names_to = "metric", values_to = "value") |> 
  mutate(value = ifelse(is.na(value), 0, value)) |> 
  pivot_wider(names_from=sensor, values_from = value, values_fill=0) |> 
  mutate(difference = PC - ARU) |> 
  rename(species_code = spp) |> 
  left_join(spp) |> 
  mutate(difference = ifelse(difference > 100, 100, difference))

write.csv(wide, file.path(root, "results", "OccurrenceSummary.csv"), row.names=FALSE)
  
#2. Plot ----
ggplot(wide) + 
  geom_jitter(aes(y=difference, x=Family, colour=Family)) +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  facet_wrap(~metric, scales="free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

#3. Test ----
  
  
  
  
  
  
  
  


  