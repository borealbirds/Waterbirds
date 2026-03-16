# ---
# title: Waterbird models - Detectability - Upland Truncation
# author: Elly Knight
# created: March 8, 2026
# ---

#NOTES################################

#PURPOSE: This script runs models to filter out upland surveys from the BAM dataset for waterbird modelling.

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
library(dismo) #BRTs
library(mgcv) #GAMs
library(pROC) #AUC for model comparison

#2. Set root path----
root <- "G:/Shared drives/BAM_WaterbirdModels"
root_data <- "G:/Shared drives/BAM_AvianData"

#3. Load waterbird data object ----
load(file.path(root, "data", "reporting_mar26", "02_WaterBirdData_Unknowns.Rdata"))

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
spp <- spp_wt |> 
  inner_join(spp_avi) |> 
  inner_join(data.frame(species_code = colnames(use))) |> 
  dplyr::select(all_of(c(colnames(spp_wt), "Family")))

#GET FLYOVER DATA ###############

#This section could probably be less jenky....

#1. Get the WildTrax data ----
load(file.path(root_data, "BAMDataset", "WildTrax", "2026-03-09", "01_wildtrax_raw_2026-03-09.Rdata"))
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
#Just waterbirds from the same surveys
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
          mutate(flyover = 1)) |> 
  left_join(spp |> 
              dplyr::select(species_code, Family))

#9. Save out for GEE ----
write.csv(all, file.path(root, "data", "WaterbirdFlyovers.csv"), row.names = FALSE)

#MODEL PREP ############

#1. Get the covariates ----
cov_pt <- read.csv(file.path(root, "data", "covariates", "WaterbirdFlyovers_water_point.csv")) |> 
  dplyr::select(-c(system.index, .geo, count)) |> 
  rename(dry_pt = b1,
         class_pt = b1_1,
         occur_pt = occurrence,
         season_pt = seasonality,
         recur_pt = recurrence) |> 
  dplyr::select(-c(b1_2, change_abs, change_norm, transition, max_extent))

cov_2km <- read.csv(file.path(root, "data", "covariates", "WaterbirdFlyovers_water_2km.csv")) |> 
  dplyr::select(-c(system.index, .geo, count)) |> 
  rename(dry_2 = b1,
         class_2 = b1_1,
         occur_2 = occurrence,
         season_2 = seasonality,
         recur_2 = recurrence) |> 
  dplyr::select(-c(b1_2, change_abs, change_norm, transition, max_extent))

cov_20km <- read.csv(file.path(root, "data", "covariates", "WaterbirdFlyovers_water_20km.csv")) |> 
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
  inner_join(cov_2km) |> 
  inner_join(cov_20km) |> 
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

#histograms
ggplot(all_wide) +
  geom_histogram(aes(x=val)) +
  facet_wrap(~cov, scales="free")

#4. Split into train and test ----
train <- all_cov |>
  sample_frac(0.9)
  
test <- all_cov |> 
  anti_join(train)

#MODEL - GLMS ###################

#Let's use:
#Linear: occur_2, recur_2,  season_2
#polynomial: dry_20, season_20, recur_20, occur_20

#1. Figure out polynomials ----

#season_20
m_season20_1 <- glm(flyover ~ season_20, data=all_cov, family = "binomial")
m_season20_2 <- glm(flyover ~ poly(season_20,2), data=all_cov, family = "binomial")
AIC(m_season20_1, m_season20_2) #second order

#recur_20
m_recur20_1 <- glm(flyover ~ recur_20, data=all_cov, family = "binomial")
m_recur20_2 <- glm(flyover ~ poly(recur_20,2), data=all_cov, family = "binomial")
AIC(m_recur20_1, m_recur20_2) #second order

m_occur20_1 <- glm(flyover ~ occur_20, data=all_cov, family = "binomial")
m_occur20_2 <- glm(flyover ~ poly(occur_20,2), data=all_cov, family = "binomial")
AIC(m_occur20_1, m_occur20_2) #second order

#2. Model no RE ----
glm1 <- glm(flyover ~ poly(dry_20, 2) + poly(occur_20, 2) + poly(recur_20, 2) + poly(season_20, 2) + occur_2 + recur_2 + season_2 + factor(class_pt),
              data=train,
              family="binomial",
              na.action = "na.fail")

d1 <- dredge(glm1)

glm2 <- glm(flyover ~ poly(dry_20, 2) + poly(occur_20, 2) + poly(recur_20, 2) + poly(season_20, 2) + season_2 + factor(class_pt),
            data=train,
            family="binomial",
            na.action = "na.fail")
summary(glm2)

#3. Model with RE for family ----
#I think this might want a random slope too
glm3 <- glmer(flyover ~ poly(dry_20, 2) + poly(occur_20, 2) + poly(recur_20, 2) + poly(season_20, 2) + occur_2 + recur_2 + season_2 + factor(class_pt) + (1|Family),
              data=train,
              family="binomial",
              na.action = "na.fail")
summary(glm3)

#4. Make predictions ----
train_pred <- data.frame(pred3 = predict(glm3, re.form = NA, type="response"),
                       pred3re = predict(glm3, type="response"),
                       pred2 = predict(glm2, type="response"),
                       pred = predict(glm1, type="response")) |> 
  cbind(train)

test_pred <- data.frame(pred3 = predict(glm3, newdata=test, re.form = NA, type="response"),
                        pred3re = predict(glm3, newdata=test, type="response"),
                        pred2 = predict(glm2, newdata=test, type="response"),
                        pred = predict(glm1, newdata=test, type="response")) |> 
  cbind(test)

#MODEL - MACHINE LEARNING #############

#1. Subset the data ----
all_brt <- train |> 
  dplyr::select(flyover, class_pt, dry_2:season_20, Family) |> 
  dplyr::select(-class_2, -class_20) |> 
  mutate(Family = as.factor(Family))

#2. Model ----
set.seed(1234)
m.i <- dismo::gbm.step(data=all_brt,
                       gbm.x=c(2:ncol(all_brt)),
                       gbm.y="flyover",
                       tree.complexity = 3,
                       learning.rate = 0.01,
                       family="bernoulli")

#3. Make predictions ----
train_pred$predbrt <- dismo::predict(m.i, x = train, type="response")
test_pred$predbrt <- dismo::predict(m.i, newdata = test, type="response")

#4. Look at some things ----
summary(m.i)
gbm.plot(m.i, smooth=TRUE)
m.int <- gbm.interactions(m.i)
par(mfrow=c(1,1))

#MODEL - GAM ###############

#1. Subset the data ----
all_gam <- all_cov |> 
  dplyr::select(flyover, dry_20, occur_20, recur_20, season_20, season_2, class_pt, Family) |> 
  mutate(Family = as.factor(Family),
         flyover = as.factor(flyover))

#2. Model ----
gam1 <- bam(flyover ~ s(dry_20, k=3) + s(occur_20, k=7) + s(recur_20, k=7) + s(season_20, k=6) + s(season_2, k=4) +s(class_pt, bs="re") + s(Family, bs="re"),
            data = all_gam,
            family = binomial(),
            gamma = 0.5)

#3. Make predictions ----
train_pred$predgam <- predict(gam1, type="response", newdata=train)
test_pred$predgam <- predict(gam1, type="response", newdata=test)

#4. Look at some things ----
summary(gam1)
gam.check(gam1)
plot.gam(gam1, pages = 1, se=FALSE)
plot.gam(gam1, select=5, rug=TRUE, ylim=c(-3, 3))

#MODEL SELECTION ##################

#1. Logloss function ----
logloss <- function(y, p){
  eps <- 1e-15
  p <- pmin(pmax(p, eps),  1 - eps)
  -mean(y * log(p) + (1-y) * log(1-p))
}

#2. Results ----
perf <- data.frame(model = rep(c("GLM", "GLMM_0", "GLMM_RE", "BRT", "GAM"), 2),
                   dataset = c(rep("train", 5), rep("test", 5)),
                   auc = c(roc(train_pred$flyover, train_pred$pred2)$auc,
                           roc(train_pred$flyover, train_pred$pred3)$auc,
                           roc(train_pred$flyover, train_pred$pred3re)$auc,
                           roc(train_pred$flyover, train_pred$predbrt)$auc,
                           roc(train_pred$flyover, train_pred$predgam)$auc,
                           roc(test_pred$flyover, test_pred$pred2)$auc,
                           roc(test_pred$flyover, test_pred$pred3)$auc,
                           roc(test_pred$flyover, test_pred$pred3re)$auc,
                           roc(test_pred$flyover, test_pred$predbrt)$auc,
                           roc(test_pred$flyover, test_pred$predgam)$auc),
                   ll = c(logloss(train_pred$flyover, train_pred$pred2),
                          logloss(train_pred$flyover, train_pred$pred3),
                          logloss(train_pred$flyover, train_pred$pred3re),
                          logloss(train_pred$flyover, train_pred$predbrt),
                          logloss(train_pred$flyover, train_pred$predgam),
                          logloss(test_pred$flyover, test_pred$pred2),
                          logloss(test_pred$flyover, test_pred$pred3),
                          logloss(test_pred$flyover, test_pred$pred3re),
                          logloss(test_pred$flyover, test_pred$predbrt),
                          logloss(test_pred$flyover, test_pred$predgam)))
perf

#3. Put together ----
all_pred <- test_pred |> 
  mutate(dataset="test") |> 
  rbind(train_pred |> 
          mutate(dataset="train"))

#3. Plot predictions ----
wide_pred <- all_pred |> 
  dplyr::select(flyover, dataset, pred3, pred3re, pred2, pred, predbrt, predgam) |> 
  pivot_longer(-c(flyover, dataset), names_to="model", values_to="prediction")

ggplot(wide_pred) +
  geom_jitter(aes(x=prediction, y=flyover, colour=dataset)) +
  geom_smooth(aes(x=prediction, y=flyover)) +
  facet_wrap(~model, ncol=2, scales="free")

#BRT seems better but concerned about overfitting.... let's try the GAM

#TRUNCATION VALUE ################

#1. Inspect ----
hist(all_pred$predbrt)

#2. Use a Gaussian mixture model to find a truncation ----
set.seed(1234)
gmm <- flexmix(predbrt ~ 1 + Family, data = all_pred, k=2, model = FLXMRglm(family = "gaussian"))
summary(gmm)

#3. Get the membership probabilities ----
all_pred$gmm_p <- posterior(gmm)[,1]
all_pred$gmm_fo <- ifelse(all_pred$gmm_p > 0.90, 1, 0)
all_pred$brt_fo <- ifelse(all_pred$predbrt > 0.2, 1, 0) #try the BRT with a cutoff

#4. Look at some things ----
table(all_pred$flyover, all_pred$gmm_fo, all_pred$dataset)  #this is garbage
table(all_pred$flyover, all_pred$brt_fo, all_pred$dataset)
ggplot(all_pred) +
  geom_jitter(aes(x=predgam, y=gmm_p, pch=factor(flyover), colour = factor(gmm_fo)), size=3) +
  facet_wrap(~flyover)

#5. Let's pick a truncation value for the BRT ----
thresh <- seq(0.1, 0.80, 0.1)
perf_thresh <- data.frame()
for(i in 1:length(thresh)){
  
  pred.i <- all_pred |> 
    mutate(fo = ifelse(predbrt > thresh[i], 1, 0)) |> 
    group_by(Family, dataset, flyover, fo) |> 
    summarize(n = n()) |> 
    ungroup() |> 
    mutate(metric = case_when(flyover==0 & fo==0 ~ "tn",
                              flyover==0 & fo==1 ~ "fp",
                              flyover==1 & fo==0 ~ "fn",
                              flyover==1 & fo==1 ~ "tp")) |> 
    dplyr::select(-flyover, -fo) |> 
    pivot_wider(names_from = metric, values_from=n)
  
  perf.i <- pred.i |> 
    group_by(Family, dataset) |> 
    summarize(recall=tp/(tp+fn),
              precision=tp/(tp+fp)) |> 
    ungroup() |> 
    mutate(thresh = thresh[i],
           fscore = (recall*precision)/(recall + precision))
  
  perf_thresh <- rbind(perf_thresh, perf.i)
  
}

ggplot(perf_thresh) + 
  geom_line(aes(x=thresh, y=recall, colour=dataset)) +
  facet_wrap(~Family)

ggplot(perf_thresh) + 
  geom_line(aes(x=thresh, y=precision, colour=dataset)) +
  facet_wrap(~Family)

ggplot(perf_thresh) + 
  geom_line(aes(x=thresh, y=fscore, colour=dataset)) +
  facet_wrap(~Family)

#CONCLUDE ##############

#1. Save a bunch of things ----
brt <- m.i
save(spp, perf_thresh, perf, all_pred, brt, gam1, glm1, glm2, glm3,
     file = file.path(root, "results", "UplandProbabilityOutput.Rdata"))

#2. Save smaller objects for reporting ----
write.csv(all_pred, file.path(root, "data", "reporting_mar26", "03_FlyoverProbability.csv"), row.names = FALSE)
write.csv(perf, file.path(root, "results", "FlyoverModelPerformance.csv"), row.names = FALSE)