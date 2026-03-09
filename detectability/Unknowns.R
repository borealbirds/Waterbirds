# ---
# title: Waterbird models - Detectability - Upland Truncation
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

#3. Load data object ----
load(file.path(root_data, "BAMDataset", "WildTrax", "2026-03-02", "01_wildtrax_raw_2026-03-02.Rdata"))

#4. Authenticate WildTrax ----
source("WTlogin.R")
wt_auth()

#5. Set unknown IDs ----
unid <- data.frame(species_code = c("UNCO", "UDOW", "UGOL", "UMAB", "UNGA", "UNLO", "UNSC", "UNWT", "UNYE", "UPHL", "USCT", "USLD", "UTEA", "UWCG"),
                   family = c("Rallidae", "Scolopacidae", "Anatidae", "Anatidae", "Rallidae", "Gaviidae", "Anatidae", "Anatidae", "Scolopacidae", "Scolopacidae", "Anatidae", "Scolopacidae", "Anatidae", "Podicipedidae"))

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
  dplyr::filter(species_code %in% unid)

#3. Proportions ----
#ARU
nrow(aru_un)/nrow(aru)*100
#0.096%

#PC
nrow(pc_un)/nrow(pc)*100
#0.004%

#CLASSIFIER RESULTS ###########

#1. Projects with ARU unknowns----
proj <- unique(aru_un$project_id)

#2. Go get the classifier report for the projects with unknowns ----
aru.list <- list()
for(i in 1:length(proj)){
  
  #authenticate each time because this loop takes forever
  wt_auth()
  
  aru.list[[i]] <- try(wt_download_report(project_id = proj[i], sensor_id = "ARU", report = "ai"))

  print(paste0("Finished dataset ", i, " of ", length(proj), " projects"))
  
}

#3. Collapse ----
ai <- do.call(rbind, aru.list)
save(ai, file = file.path(root, "data", "WaterbirdFlyovers_ClassifierReports.Rdata"))

#4. Wrangle to recordings of interest ----
ai_rec <- ai |> 
  dplyr::filter(recording_id %in% aru_un$recording_id,
                version=="Hawkears v1.0.8")
rm(ai, aru.list)

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

#4. Quantify ----
table(aru_ai$ai_match, aru_ai$he_report)

#number of potential classifications - 62.5%
ai_total <- aru_ai |> 
  dplyr::filter(he_report=="Yes_HE") |> 
  mutate(unid = n()) |> 
  group_by(ai_match, unid) |> 
  summarize(nmatch = n()) |> 
  ungroup() |> 
  mutate(percent = nmatch/unid*100)
ai_total

#remaining unid percentage - 0.033%
ai_total |> 
  dplyr::filter(ai_match=="No_match") |> 
  summarize(percent = nmatch/nrow(aru)*100)
#still 10 times as much as point count, but pretty negligible

#5. Look at amplitude ----
#Are they just really quiet detections that we probably wouldn't get on a point count anyway?
ggplot(aru_ai) +
  geom_histogram(aes(x=rms_peak_dbfs, fill=ai_match))

t.test(rms_peak_dbfs ~ ai_match, data=aru_ai)
