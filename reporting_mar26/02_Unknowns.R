# ---
# title: Waterbird models - Unknowns
# author: Elly Knight
# created: March 8, 2026
# ---

#NOTES################################

#PURPOSE: This script determines the proportion of unknown waterbird species identifications in acoustic recordings and explores the use of HawkEars for reduction of unknowns.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(terra) #rasters
library(sf) #shapefiles
library(wildrtrax) #species codes list
library(avilistr) #taxonomy

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"
root_data <- "G:/Shared drives/BAM_AvianData"

#3. Load data objects ----
load(file.path(root_data, "BAMDataset", "WildTrax", "2026-03-09", "01_wildtrax_raw_2026-03-09.Rdata"))
load(file.path(root, "data", "reporting_mar26", "01_WaterBirdData_Habitat.Rdata"))

#4. Authenticate WildTrax ----
source("WTlogin.R")
wt_auth()

#5. Set unknown IDs ----
unid <- data.frame(species_code = c("UNCO", "UNQK", "UGRB", "UNME", "UNGO", "UDAB", "UDIV", "UDOW", "UGOL", "UMAB", "UNDU",  "UNGA", "UNLO", "UNSC", "UNWT", "UNYE", "UPHL", "USCT", "USLD", "UTEA", "UWCG", "UNSH", "USWN"),
                   family = c("Rallidae", "Anatidae", "Podicipedidae", "Anatidae", "Anatidae", "Anatidae", "Anatidae", "Scolopacidae", "Anatidae", "Anatidae", "Anatidae", "Rallidae", "Gaviidae", "Anatidae", "Anatidae", "Scolopacidae", "Scolopacidae", "Anatidae", "Scolopacidae", "Anatidae", "Podicipedidae", "Charadriidae", "Anatidae"))

#WRANGLING###########

#1. Rbind the objects ----
aru <- do.call(rbind, aru.wt) |> 
  dplyr::select(-max_noise_type, -max_noise_volume, -max_noise_channel, -max_noise_density) |> 
  unique()
pc <- do.call(rbind, pc.wt)

#2. Get the unknowns ----
aru_un <- aru |> 
  dplyr::filter(species_code %in% unid$species_code) |> 
  rename(tag_start = detection_time) |> 
  mutate(tag_end = tag_start + tag_duration)

pc_un <- pc |> 
  dplyr::filter(species_code %in% unid$species_code)

#3. Proportions ----
#ARU
nrow(aru_un)/nrow(aru)*100
#0.36%

#PC
nrow(pc_un)/nrow(pc)*100
#0.014%

#4. Save for reporting ----
un <- aru_un |> 
  dplyr::select(organization, location_id, latitude, longitude, species_code, species_common_name, rms_peak_dbfs) |> 
  mutate(sensor="ARU") |> 
  rbind(pc_un |> 
          dplyr::select(organization, location_id, latitude, longitude, species_code, species_common_name) |> 
          mutate(rms_peak_dbfs = NA,
                 sensor="PointCount"))
write.csv(un, file.path(root, "data", "WaterbirdUnknowns.csv"), row.names=FALSE)

#CLASSIFIER RESULTS ###########

#This takes a long time, only run it once

# #1. Projects with ARU unknowns----
# proj <- unique(aru_un$project_id)
# 
# #2. Go get the classifier report for the projects with unknowns ----
# aru.list <- list()
# for(i in 1:length(proj)){
# 
#   #authenticate each time because this loop takes forever
#   wt_auth()
# 
#   aru.list[[i]] <- try(wt_download_report(project_id = proj[i], sensor_id = "ARU", report = "ai"))
# 
#   print(paste0("Finished dataset ", i, " of ", length(proj), " projects"))
# 
# }
# 
# #3. Collapse ----
# ai <- do.call(rbind, aru.list)
# save(ai, file = file.path(root, "data", "reporting_mar26", "WaterbirdFlyovers_ClassifierReports.Rdata"))
# load(file.path(root, "data", "reporting_mar26", "WaterbirdFlyovers_ClassifierReports.Rdata"))

#4. Wrangle to recordings of interest ----
ai_rec <- ai |> 
  dplyr::filter(recording_id %in% aru_un$recording_id,
                version=="Hawkears v1.0.8")
rm(ai, aru.list)

write.csv(ai_rec, file.path(root, "data", "reporting_mar26", "WaterbirdFlyovers_ClassifierRecordings.csv"), row.names=FALSE)
ai_rec <- read.csv(file.path(root, "data", "reporting_mar26", "WaterbirdFlyovers_ClassifierRecordings.csv"))

#RECOMMENDATIONS ##########

#1. Filter to detections in the same time window of the same family ----
ai_un <- ai_rec |> 
  rename(species_ai = species_code) |> 
  dplyr::filter(species_ai %in% spp$species_code,
                is_species_allowed_in_project==TRUE) |> 
  left_join(aru_un |> 
              dplyr::select(recording_id, tag_id, tag_start, tag_end, species_code)) |> 
  dplyr::filter((detection_time > (tag_start-3) & detection_time < tag_end) |
                  (detection_time > tag_start) & detection_time < (tag_end+3)) |> 
  unique() |> 
  left_join(spp |> 
              rename(species_ai = species_code,
                     family_ai = Family) |> 
              dplyr::select(species_ai, family_ai)) |> 
  left_join(unid) |> 
  dplyr::filter(family_ai==family) |> 
  group_by(recording_id, tag_id, tag_start, tag_end, species_code, species_ai, family) |> 
  summarize(score = max(confidence)) |> 
  ungroup()

#2. Pick the highest scoring species for tags with multiple options in the family ----
ai_max <- ai_un |> 
  group_by(recording_id, tag_id, tag_start, tag_end, species_code, family) |> 
  mutate(maxscore = max(score),
         ai_tags = n()) |> 
  arrange(tag_id, -score) |> 
  ungroup() |> 
  mutate(scorediff = ifelse(ai_tags > 1 & score==maxscore, score - lead(score), NA)) |> 
  dplyr::filter(score==maxscore) |> 
  dplyr::select(-maxscore)

#3. Join back to the rest of the unknowns ----
aru_ai <- aru_un |> 
  left_join(ai_max |> 
              dplyr::select(tag_id, species_ai, score, ai_tags, scorediff)) |> 
  mutate(ai_match = ifelse(is.na(species_ai), "No_match", "AI_match"),
         he_report = ifelse(recording_id %in% ai_rec$recording_id, "Yes_HE", "No_HE"))

#4. More cleaning ----
ai_clean <- aru_ai |> 
  mutate(species_replace = case_when(species_code=="UDAB" & species_ai %in% c("GADW", "MALL") ~ species_ai,
                                     species_code=="UDIV" & species_ai %in% c("BUFF", "COGO", "RNDU") ~ species_ai,
                                     species_code=="UGOL" & species_ai %in% c("COGO", "BAGO") ~ species_ai,
                                     species_code=="UGRB" & species_ai %in% c("PBGR") ~ species_ai,
                                     species_code=="UMAB" & species_ai %in% c("Mall", "ABDU") ~ species_ai,
                                     species_code=="UNDU" & !species_ai %in% c("CANG", "BRAN", "GWFG", "SNGO", "TRUS", "TUSW") ~ species_ai,
                                     species_code=="UNGO" & species_ai %in% c("CANG", "BRAN", "GWFG", "SNGO", "TRUS", "TUSW") ~ species_ai,
                                     species_code=="UNLO" ~ species_ai,
                                     species_code=="UNSC" ~ NA,
                                     species_code=="UNYE" & species_ai %in% c("GRYE", "LEYE") ~ species_ai,
                                     species_code=="USWN" & species_ai %in% c("TRUS", "TUSW") ~ species_ai,
                                     species_code=="UTEA" & species_ai %in% c("GWTE", "BLTE", "CITE") ~ species_ai,
                                     species_code=="UWCG" & species_ai %in% c("WEGR", "CLGR") ~ species_ai)) |>
  mutate(ai_clean = ifelse(is.na(species_replace), "No_match", "Clean_match"))

write.csv(ai_clean, file.path(root, "data", "reporting_mar26", "WaterbirdUnknowns_HawkEarsReplacement.csv"), row.names=FALSE)

#5. Data presents ----
table(ai_clean$ai_match, ai_clean$he_report)
table(ai_clean$ai_match, ai_clean$ai_clean)
table(ai_clean$species_code, ai_clean$ai_clean)
table(ai_clean$species_ai, ai_clean$species_code)

ggplot(ai_clean) +
  geom_histogram(aes(x=score, fill=ai_clean))

ggplot(ai_clean |> dplyr::filter(ai_clean=="Clean_match")) +
  geom_histogram(aes(x=score)) +
  facet_wrap(~species_code, scales="free")

ggplot(ai_clean |> dplyr::filter(ai_clean=="Clean_match")) +
  geom_point(aes(x=rms_peak_dbfs, y=score)) +
  geom_smooth(aes(x=rms_peak_dbfs, y=score))

#Number of recordings with HE
ai_clean |> 
  group_by(he_report) |> 
  summarize(n =n()) |> 
  ungroup() |> 
  mutate(percent = n/nrow(ai_clean)*100)


#number of potential classifications - 62.5%
ai_total <- ai_clean |> 
  dplyr::filter(he_report=="Yes_HE") |> 
  mutate(unid = n()) |> 
  group_by(ai_clean, unid) |> 
  summarize(nmatch = n()) |> 
  ungroup() |> 
  mutate(percent = nmatch/unid*100)
ai_total

#remaining unid percentage - 0.033%
ai_total |> 
  dplyr::filter(ai_clean=="No_match") |> 
  summarize(percent = nmatch/nrow(aru)*100)
#still 10 times as much as point count, but pretty negligible

#5. Look at amplitude ----
#Are they just really quiet detections that we probably wouldn't get on a point count anyway?
ggplot(ai_clean) +
  geom_histogram(aes(x=rms_peak_dbfs, fill=ai_clean))

t.test(rms_peak_dbfs ~ ai_clean, data=ai_clean)
t.test(score ~ ai_clean, data=ai_clean)

#DATASET ADJUSTMENT #############

#1. Filter to just replacement ones ----
ai_replace <- ai_clean |> 
  dplyr::filter(ai_clean=="Clean_match") |> 
  rename(survey_id = task_id) |> 
  dplyr::select(survey_id, species_replace) |> 
  group_by(survey_id, species_replace) |> 
  summarize(add = n()) |> 
  ungroup()

#2. Get the surveys that need adjustment ----
adjust <- dat |> 
  mutate(survey_id = as.numeric(survey_id)) |> 
  dplyr::filter(survey_id %in% ai_replace$survey_id) |> 
  pivot_longer(LBDO:CITE, names_to = "species_replace", values_to="count")

#3. Put together and make it long again ----
adjust_ai <- adjust |> 
  left_join(ai_replace) |> 
  mutate(add = ifelse(is.na(add), 0, add),
         newcount = count + add)  |> 
  dplyr::select(-count, -add) |> 
  pivot_wider(names_from = species_replace, values_from=newcount)

#4. Add back to dataset ----
clean <- dat |> 
  dplyr::filter(!survey_id %in% adjust_ai$survey_id) |> 
  rbind(adjust_ai)

#5. Save ----
save(clean, cov, spp, file = file.path(root, "data", "reporting_mar26", "02_WaterBirdData_Unknowns.Rdata"))
