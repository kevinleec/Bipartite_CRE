##### Load Libraries & Data #####
library(gdata)
library(stringr)
library(dplyr)
library(data.table)
library(xgboost)
library(ggplot2)

# load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/HyADSmat.Rda")
# load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/out.zip_pp.rda")
# load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/facilities_for_analysis.Rda")
# load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/ZipcodeData.Rda")

load("~/OneDrive - Harvard University/Research/Bipartite CRE/Simulations/sim_data.RData")
load("~/OneDrive - Harvard University/Research/Bipartite CRE/Simulations/HyADSmat.Rda")
load("~/OneDrive - Harvard University/Research/Bipartite CRE/Simulations/out.zip_pp.rda")
load("~/OneDrive - Harvard University/Research/Bipartite CRE/Simulations/facilities_for_analysis.RDa")

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

## set power plant and zipcode covariates to be included
covs <- c("logPop", "PctUrban", "PctHisp", "PctHighSchool", "logMedianHHInc", "PctPoor", "PctOccupied", "PctMovedIn5", "logPopPerSQM", #census
          "smokerate", #smoking
          "temp", "rh", #weather
          "totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2", #power plant characteristics
          "mean_age", "Female_rate", "White_rate", "Black_rate", #Medicare characteristics
          "HyADS") #HyADS
zip.covs <- c("logPop", "PctUrban", "PctHisp", "PctHighSchool", "logMedianHHInc", "PctPoor", "PctOccupied", "PctMovedIn5", "logPopPerSQM", #census
              "smokerate", #smoking
              "temp", "rh", #weather
              "mean_age", "Female_rate", "White_rate", "Black_rate") #Medicare characteristics
pp.covs = c("totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2")  # - these are power plant covariates that are very predictive in PS model, but likely not confounders

downweight_ratio <- 0.6786576 ## mean of ratio of 2nd key-associated and key-associated HyADS for each zip code

################################################################
## Simulation with Key-Associated Plants
################################################################

##### Data Preparation and Cleaning #####
hyads_var <- as.data.frame(sumMAP.2005.fac) ## hyads matrix
zip_var <- data ## zipcode covariates
facility_var_2005 <- dat_facility_2005 ## power plant covariates
rm("sumMAP.2005.fac")
rm("data")
rm("dat_facility_2005")


zips <- data.frame("ZIP" = zip_var$ZIP) ## list of zipcodes
plants <- data.frame("FacID" = facility_var_2005$FacID) ## list of power plant IDs

## transform covariates
zip_var$logPop <- log(zip_var$TotPop)
zip_var$logMedianHHInc <- log(zip_var$MedianHHInc)
zip_var$logPopPerSQM <- log(zip_var$PopPerSQM)
facility_var_2005$logHeatInput <- log(facility_var_2005$totHeatInput)
facility_var_2005$logOpTime <- log(facility_var_2005$totOpTime)

zip_covariates <- select(zip_var, c("ZIP", zip.covs)) ## select zipcode covariates to be used
plant_covariates <- select(facility_var_2005, c("FacID", pp.covs)) ## select power plant covariates to be used


##### Without filtering for complete covariate data #####
## connect key-associated and 2nd key-associated plant to each outcome unit
dat <- as.data.frame(t(sapply(hyads_var, function(x) head(row.names(hyads_var)[order(x, decreasing = TRUE)], 2))))
setDT(dat, keep.rownames = "ZIP")
colnames(dat) <- c("ZIP", "Key", "SecondKey")
dat <- as.data.frame(dat)
#dat$Z <- as.integer(dat$Z)
dat$Key <- as.character(dat$Key)
dat$SecondKey <- as.character(dat$SecondKey)
dat$Key_HyADS <- hyads_var[as.matrix(cbind(dat$Key, dat$ZIP))] ## hyads for key-associated plants
dat$SecondKey_HyADS <- hyads_var[as.matrix(cbind(dat$SecondKey, dat$ZIP))]

## aggregate mean of covariates for each power plant for its key-associated zipcodes
zip_agg <- dat
zip_agg <- zip_agg %>% ## merge key-associated with zip covariates
  full_join(zip_covariates, by = "ZIP")
zip_agg[zip_agg == -Inf] <- NA ## convert -Inf to NA

## create summary of zipcodes key-associated and second key-associated with plant
key_covariates <- aggregate(zip_agg[,-c(1:5)], list(zip_agg$Key), mean, na.rm=T)
colnames(key_covariates) <- c("Key", paste("Key", colnames(key_covariates)[-1], sep = "_"))
secondkey_covariates <- aggregate(zip_agg[,-c(1:5)], list(zip_agg$SecondKey), mean, na.rm=T)
colnames(secondkey_covariates) <- c("SecondKey", paste("SecondKey",colnames(secondkey_covariates)[-1], sep = "_"))
secondkey_covariates[,-1] <- secondkey_covariates[,-1] * downweight_ratio ## downweight second key-associated

plant_dat <- key_covariates %>% ## plant, key-associated covariates, second key-associated covariates, plant covariates
  full_join(secondkey_covariates, by = c("Key"="SecondKey")) %>%
  left_join(plant_covariates, by = c("Key"="FacID"))





##### Filter by complete covariate data #####
complete_second <- plant_dat[complete.cases(plant_dat[,18:33]),]
complete_key <- plant_dat[complete.cases(plant_dat[,2:17]),]
complete_plants <- plant_dat[complete.cases(plant_dat[,34:39]),]
complete_all <- plant_dat[complete.cases(plant_dat),]
plants_complete <- complete_all$Key


hyads_mat_complete <- hyads_var[plants_complete,] ## hyads matrix based on plants with complete data
## connect key-associated and 2nd key-associated plant to each outcome unit
dat_complete <- as.data.frame(t(sapply(hyads_mat_complete, function(x) head(row.names(hyads_mat_complete)[order(x, decreasing = TRUE)], 2))))
setDT(dat_complete, keep.rownames = "ZIP")
colnames(dat_complete) <- c("ZIP", "Key", "SecondKey")
dat_complete <- as.data.frame(dat_complete)
#dat$Z <- as.integer(dat$Z)
dat_complete$Key <- as.character(dat_complete$Key)
dat_complete$SecondKey <- as.character(dat_complete$SecondKey)
dat_complete$Key_HyADS <- hyads_mat_complete[as.matrix(cbind(dat_complete$Key, dat_complete$ZIP))] ## hyads for key-associated plants
dat_complete$SecondKey_HyADS <- hyads_mat_complete[as.matrix(cbind(dat_complete$SecondKey, dat_complete$ZIP))]

changes <- anti_join(dat, dat_complete, by=c("ZIP", "Key", "SecondKey")) ## look at changes

## generate treatments
treatments <- left_join(complete_all, facility_var_2005[,c(1,7)], by = c("Key"="FacID"))
# treat.mod <- glm(ScrubbedFacility ~ . - Key, test, family = binomial())
# test$prop <- predict(treat.mod, newdata=test, type = "response")

treatments <- treatments %>% mutate(x = -100 - .1*Key_logPop + 1*Key_PctUrban + 5*Key_PctHisp - 1*Key_PctHighSchool - 1*Key_logMedianHHInc - 10*Key_PctPoor - 5*Key_PctOccupied + 1*Key_PctMovedIn5 - 2*Key_logPopPerSQM - 2*Key_smokerate + .1*Key_temp + 500*Key_rh + 1*Key_mean_age - 10*Key_Female_rate - 8*Key_White_rate - 9*Key_Black_rate +
  1*SecondKey_logPop + 5*SecondKey_PctUrban -10*SecondKey_PctHisp + 1*SecondKey_PctHighSchool + 1*SecondKey_logMedianHHInc + 20*SecondKey_PctPoor - 2*SecondKey_PctOccupied - 2*SecondKey_PctMovedIn5 - 1*SecondKey_logPopPerSQM - 2*SecondKey_smokerate + .1*SecondKey_temp - 1000*SecondKey_rh + .5*SecondKey_mean_age - 10*SecondKey_Female_rate + 1*SecondKey_White_rate + 1*SecondKey_Black_rate +
  .5*totNumNOxControls - 1*logHeatInput - 1*logOpTime + 5*pctCapacity + 1*pctS_n_CR - 1*Phase2)

treatments$prop <- expit(treatments$x)
treatments$Z <- rbinom(nrow(treatments),1,treatments$prop)
complete_all <- full_join(complete_all, treatments[,c(1,43)], by = "Key")

## merge zip-level data with key-associated plant covariates and covariate summaries
dat_complete <- dat_complete %>%
left_join(plant_covariates, by = c("Key"="FacID")) %>%
left_join(key_covariates, by = "Key") %>%
left_join(secondkey_covariates, by = "SecondKey")

## Get key and upwind treatments using plant treatment mapping
treat_mapping <- setNames(complete_all$Z, complete_all$Key)
dat_complete$Z <- as.integer(treat_mapping[dat_complete$Key]) ## based on key-associated plant: Z=1 if treated, Z=0 if not
dat_complete$G <- as.integer(treat_mapping[dat_complete$SecondKey])

## remove treatment NAs
#dat <- subset(dat, !is.na(Z) & !is.na(G)) ## unnecessary for filtered data

## generate outcome
dat_complete <- left_join(dat_complete, zip_covariates, by = "ZIP")
dat_complete <- dat_complete[complete.cases(dat_complete),]
dat_complete <- dat_complete %>% mutate(lambda = -100*Z - 67*G + 100*logPop + 80*PctUrban + 50*PctHisp - 80*PctHighSchool
  - 5*logMedianHHInc + 100*PctPoor - 2*PctOccupied + 1*PctMovedIn5
  + 5*logPopPerSQM + 200*smokerate + .1*temp - 1000*rh + 1*mean_age
  - 30*Female_rate - 70*White_rate + 150*Black_rate)

# dat_complete <- dat_complete %>% mutate(loglambda = -4.6*Z - 3*G + 4.6*logPop + 3*PctUrban + 3.9*PctHisp - 4.4*PctHighSchool
#                                         - 1.6*logMedianHHInc + 4.6*PctPoor - 0.7*PctOccupied + 0.1*PctMovedIn5
#                                         + 1.6*logPopPerSQM + 5.3*smokerate - .02*temp - 7*rh + 0.1*mean_age
#                                         - 3.4*Female_rate - 4.6*White_rate + 5*Black_rate)
# dat_complete$lambda <- exp(dat_complete$loglambda)

dat_complete$Y <- rpois(nrow(dat_complete), dat_complete$lambda)
colnames(dat_complete)
dat_complete <- dat_complete[,-c(12:43)] ## remove power plant-level summarized covariates


## histogram of outcome by group; means of each group (Z=z,G=g)
ggplot(dat_complete, aes(x=Y)) +
  geom_histogram(color="#e9ecef", aes(fill=as.factor(Z)), position="identity", bins=30, alpha=0.6) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="Z")

mean(dat_complete$Y[dat_complete$Z==1 & dat_complete$G==1])
mean(dat_complete$Y[dat_complete$Z==1 & dat_complete$G==0])
mean(dat_complete$Y[dat_complete$Z==0 & dat_complete$G==1])
mean(dat_complete$Y[dat_complete$Z==0 & dat_complete$G==0])




##### Propensity score estimation #####

Xzip.formula <- paste(zip.covs, collapse=" + ")
formulaZ <- as.formula(paste('Z~', Xzip.formula))
formulaG <- as.formula(paste('G~Z+', Xzip.formula))

## Estimate the key associated propensity score - Z | X_out
psmodZ <- glm(formulaZ, dat=subset(dat_complete, Z %in% c(0,1)), family=binomial(link="logit"))
summary(psmodZ)
dat_complete <- as.data.table(dat_complete)
dat_complete[Z%in%c(0,1), psZ := psmodZ$fitted.values]
#dat_complete[Z%in%c(0,1), psZ_logit := psmodZ$linear.predictors]

# ## xgboost
# train_prop <- 0.75
# set.seed(123)
# boostZ_covariates <- c("Key_HyADS", "SecondKey_HyADS", zip.covs, "Z")
# datZ <- select(dat_complete, all_of(boostZ_covariates))
# 
# train_indZ <- sample(seq_len(nrow(datZ)), size = floor(train_prop * nrow(datZ)))
# df_trainZ <- as.data.frame(datZ[train_indZ, ])
# df_testZ <- as.data.frame(datZ[-train_indZ, ])
# 
# train_propZ <- df_trainZ['Z']
# train_XZ <- df_trainZ[-grep('Z', colnames(df_trainZ))]
# test_propZ <- df_testZ['Z']
# test_XZ <- df_testZ[-grep('Z', colnames(df_testZ))]
# 
# xgbcv <- xgb.cv(data = as.matrix(train_XZ), label = t(t(train_propZ)), 
#                 max_depth = 6,
#                 eta = 0.05,
#                 gamma = 5,
#                 lambda = 2,
#                 subsample = 0.5,
#                 nrounds = 100,
#                 objective = "binary:logistic",
#                 verbose = 1,
#                 nfold = 5,
#                 print_every_n = 10,
#                 early_stopping_rounds = 20,
#                 maximize = F,
#                 eval_metric = "error")
# 
# xgbModZ <- xgboost(data = as.matrix(train_XZ), label = t(t(train_propZ)), 
#                    max_depth = 6,
#                    eta = 0.05,
#                    gamma = 5,
#                    lambda = 2,
#                    subsample = 0.5,
#                    nrounds = 50,
#                    objective = "binary:logistic",
#                    verbose = 1)
# 
# prop_predZ <- predict(xgbModZ, as.matrix(test_XZ))
# 
# errZ <- mean(as.numeric(prop_predZ > 0.5) != test_propZ)
# print(paste("test-error=", errZ))
# 
# datZ <- as.data.frame(datZ)
# dat_complete$xg_psZ <- predict(xgbModZ, as.matrix(datZ[-grep('Z', colnames(datZ))]))






## Estimate upwind propensity score - G | Z, X_out
psmodG <- glm(formulaG, dat=subset(dat_complete, G %in% c(0,1)), family=binomial(link="logit"))
summary(psmodG)
dat_complete[G%in%c(0,1), psG := psmodG$fitted.values]
#dat_complete[G%in%c(0,1), psG_logit := psmodG$linear.predictors]

# ## xgboost
# boostG_covariates <- c("Key_HyADS", "SecondKey_HyADS", zip.covs, "Z", "G")
# datG <- select(dat_complete, all_of(boostG_covariates))
# 
# train_indG <- sample(seq_len(nrow(datG)), size = floor(train_prop * nrow(datG)))
# df_trainG <- as.data.frame(datG[train_indG, ])
# df_testG <- as.data.frame(datG[-train_indG, ])
# 
# train_propG <- df_trainG['G']
# train_XG <- df_trainG[-grep('G', colnames(df_trainG))]
# test_propG <- df_testG['G']
# test_XG <- df_testG[-grep('G', colnames(df_testG))]
# 
# xgbcv <- xgb.cv(data = as.matrix(train_XG), label = t(t(train_propG)), 
#                 max_depth = 6,
#                 eta = 0.05,
#                 gamma = 5,
#                 lambda = 2,
#                 subsample = 0.5,
#                 nrounds = 100,
#                 objective = "binary:logistic",
#                 verbose = 1,
#                 nfold = 5,
#                 print_every_n = 10,
#                 early_stopping_rounds = 20,
#                 maximize = F,
#                 eval_metric = "error")
# 
# xgbModG <- xgboost(data = as.matrix(train_XG), label = t(t(train_propG)), 
#                    max_depth = 6,
#                    eta = 0.05,
#                    gamma = 5,
#                    lambda = 2,
#                    subsample = 0.5,
#                    nrounds = 50,
#                    objective = "binary:logistic",
#                    verbose = 1)
# 
# prop_predG <- predict(xgbModG, as.matrix(test_XG))
# 
# errG <- mean(as.numeric(prop_predG > 0.5) != test_propG)
# print(paste("test-error=", errG))
# 
# datG <- as.data.frame(datG)
# dat_complete$xg_psG <- predict(xgbModG, as.matrix(datZ[-grep('G', colnames(datG))]))

## Joint propensity scores
dat_complete$Z1_G1 <- dat_complete$psZ * dat_complete$psG ## P(Z=1,G=1 | X_out)
dat_complete$Z1_G0 <- dat_complete$psZ * (1-dat_complete$psG) ## P(Z=1,G=0 | X_out)
dat_complete$Z0_G1 <- (1-dat_complete$psZ) * dat_complete$psG ## P(Z=0,G=1 | X_out)
dat_complete$Z0_G0 <- (1-dat_complete$psZ) * (1-dat_complete$psG) ## P(Z=0,G=0 | X_out)





##### Estimator #####



## outcome model
formulaY <- as.formula(paste('Y ~ ', Xzip.formula))

subset_Z1G1 <- subset(dat_complete, Z == 1 & G == 1)
subset_Z1G0 <- subset(dat_complete, Z == 1 & G == 0)
subset_Z0G1 <- subset(dat_complete, Z == 0 & G == 1)
subset_Z0G0 <- subset(dat_complete, Z == 0 & G == 0)
out_mod_Z1G1 <- glm(formulaY, data = subset_Z1G1, family = quasipoisson(link = "log"))
out_mod_Z1G0 <- glm(formulaY, data = subset_Z1G0, family = quasipoisson(link = "log"))
out_mod_Z0G1 <- glm(formulaY, data = subset_Z0G1, family = quasipoisson(link = "log"))
out_mod_Z0G0 <- glm(formulaY, data = subset_Z0G0, family = quasipoisson(link = "log"))


newdata_Z1G1 <- dat_complete %>% mutate(Z = 1, G = 1)
newdata_Z1G0 <- dat_complete %>% mutate(Z = 1, G = 0)
newdata_Z0G1 <- dat_complete %>% mutate(Z = 0, G = 1)
newdata_Z0G0 <- dat_complete %>% mutate(Z = 0, G = 0)

dat_complete <- dat_complete %>% 
  mutate(predY_Z1G1 = predict(out_mod_Z1G1, newdata = newdata_Z1G1, type = "response"),
         predY_Z1G0 = predict(out_mod_Z1G0, newdata = newdata_Z1G0, type = "response"),
         predY_Z0G1 = predict(out_mod_Z0G1, newdata = newdata_Z0G1, type = "response"),
         predY_Z0G0 = predict(out_mod_Z0G0, newdata = newdata_Z0G0, type = "response"))



aipw_po <- function(z,g) {
  
  aipw_df <- dat_complete %>%
    mutate(po = case_when(
      z == 1 & g == 1 ~ (1/Z1_G1)*Y + (1 - 1/Z1_G1)*predY_Z1G1,
      z == 1 & g == 0 ~ (1/Z1_G0)*Y + (1 - 1/Z1_G0)*predY_Z1G0,
      z == 0 & g == 1 ~ (1/Z0_G1)*Y + (1 - 1/Z0_G1)*predY_Z0G1,
      z == 0 & g == 0 ~ (1/Z0_G0)*Y + (1 - 1/Z0_G0)*predY_Z0G0
    )
  )
  
  return(aipw_df$po)

}
sum(aipw_po(1,0) - aipw_po(0,0) > 0)










## -----------------------------------ignore this--------------------------------- ##
## generate outcome: include treatments, zipcode-level covariates





# key.covs <- c("Key_logPop", "Key_PctUrban", "Key_PctHisp", "Key_PctHighSchool", "Key_logMedianHHInc", "Key_PctPoor", "Key_PctOccupied", "Key_PctMovedIn5", "Key_logPopPerSQM", #census
#               "Key_smokerate", #smoking
#               "Key_temp", "Key_rh", #weather
#               "Key_mean_age", "Key_Female_rate", "Key_White_rate", "Key_Black_rate") #Medicare characteristics
# secondkey.covs <- c("SecondKey_logPop", "SecondKey_PctUrban", "SecondKey_PctHisp", "SecondKey_PctHighSchool", "SecondKey_logMedianHHInc", "SecondKey_PctPoor", "SecondKey_PctOccupied", "SecondKey_PctMovedIn5", "SecondKey_logPopPerSQM", #census
#                     "SecondKey_smokerate", #smoking
#                     "SecondKey_temp", "SecondKey_rh", #weather
#                     "SecondKey_mean_age", "SecondKey_Female_rate", "SecondKey_White_rate", "SecondKey_Black_rate") #Medicare characteristics



## Joint propensity score components:
## -- Z: P( key-associated power plant treated | power plant covariates, key-associated zip code covariates )
## -- G: P( 2nd key-associated power plant treated | Z, power plant covariates, 2nd key-associated zip code covariates )

## generate intervention-level covariates and treatment
# plants$W1 <- rnorm(nrow(plants), 0, 0.2^2)
# plants$W2 <- rnorm(nrow(plants), 0, 0.2^2)
# plants$Z <- rbinom(nrow(plants), 1, expit(-0.2 + 0.3*plants$W1 - 0.15*plants$W2))

## generate upwind treatment status
# upwind_mapping <- setNames(plants$Z, plants$FacID)
# dat$G <- upwind_mapping[dat$Second_Key] ## based on second key-associated plant: G=1 if treated, G=0 if not
# dat$G <- as.integer(dat$G)

# ## generate outcome-level covariates
# dat$X1 <- rnorm(nrow(dat), 0, 0.2^2)
# dat$X2 <- rnorm(nrow(dat), 0, 0.2^2)




##### Propensity Score #####

outcome_covariates <- c("X1", "X2")
intervention_covariates <- c("W1", "W2")

Wz.formula <- paste(outcome_covariates, collapse=" + ")
formulaZ <- as.formula(paste('Z~', Wz.formula))

Wg.formula <- paste(outcome_covariates, collapse=" + ")
formulaG <- as.formula(paste('G~Z+', Wg.formula))

## Estimate the key associated propensity score - Z | W, X
psmodZ <- glm(formulaZ, dat=subset(dat, Z %in% c(0,1)), family=binomial(link="logit"))
summary(psmodZ)
dat <- as.data.table(dat)
dat[Z%in%c(0,1), psZ := psmodZ$fitted.values]
dat[Z%in%c(0,1), psZ_logit := psmodZ$linear.predictors]


## -- xgboost
train_prop <- 0.75
set.seed(123)
boostZ_covariates <- c("Key_HyADS", "Second_Key_HyADS", "W1", "W2", "Z", "X1", "X2")
datZ <- select(dat, all_of(boostZ_covariates))

train_indZ <- sample(seq_len(nrow(datZ)), size = floor(train_prop * nrow(datZ)))
df_trainZ <- as.data.frame(datZ[train_indZ, ])
df_testZ <- as.data.frame(datZ[-train_indZ, ])

train_propZ <- df_trainZ['Z']
train_XZ <- df_trainZ[-grep('Z', colnames(df_trainZ))]
test_propZ <- df_testZ['Z']
test_XZ <- df_testZ[-grep('Z', colnames(df_testZ))]

xgbModZ <- xgboost(data = as.matrix(train_XZ), label = t(t(train_propZ)), 
                  max_depth = 2,
                  eta = 1,
                  nthread = 2,
                  nrounds = 30,
                  objective = "binary:logistic",
                  verbose = 1)

prop_predZ <- predict(xgbModZ, as.matrix(test_XZ))

errZ <- mean(as.numeric(prop_predZ > 0.5) != test_propZ)
print(paste("test-error=", errZ))

datZ <- as.data.frame(datZ)
dat$xg_psZ <- predict(xgbModZ, as.matrix(datZ[-grep('Z', colnames(datZ))]))


## Estimate upwind propensity score - G | Z, W, X
psmodG <- glm(formulaG, dat=subset(dat, G %in% c(0,1)), family=binomial(link="logit"))
summary(psmodG)
dat[G%in%c(0,1), psG := psmodG$fitted.values]
dat[G%in%c(0,1), psG_logit := psmodG$linear.predictors]


## -- xgboost
train_prop <- 0.75
set.seed(123)
boostG_covariates <- c("Key_HyADS", "Second_Key_HyADS", "W1", "W2", "Z", "G", "X1", "X2")
datG <- select(dat, all_of(boostG_covariates))

train_indG <- sample(seq_len(nrow(datG)), size = floor(train_prop * nrow(datG)))
df_trainG <- as.data.frame(datG[train_indG, ])
df_testG <- as.data.frame(datG[-train_indG, ])

train_propG <- df_trainG['G']
train_XG <- df_trainG[-grep('G', colnames(df_trainG))]
test_propG <- df_testG['G']
test_XG <- df_testG[-grep('G', colnames(df_testG))]

xgbModG <- xgboost(data = as.matrix(train_XG), label = t(t(train_propG)), 
                  max_depth = 2,
                  eta = 1,
                  nthread = 2,
                  nrounds = 30,
                  objective = "binary:logistic",
                  verbose = 1)

prop_predG <- predict(xgbModG, as.matrix(test_XG))

errG <- mean(as.numeric(prop_predG > 0.5) != test_propG)
print(paste("test-error=", errG))


datG <- as.data.frame(datG)
dat$xg_psG <- predict(xgbModG, as.matrix(datG[-grep('G', colnames(datG))]))


## joint propensity score

dat$jps = dat$xg_psZ * dat$xg_psG


##### Outcome Analysis #####

## simulate outcome


## inverse joint propensity score weighted estimator

dat$Y_ipw <- (dat$Y * dat$Z * dat$G) / dat$jps




################################################################
## Long form analysis code (zip-powerplant combinations)
################################################################

##### Pseudo-simulation and cleaning #####
dat_hyads <- as.data.frame(unmatrix(hyads_var))
colnames(dat_hyads) <- "HyADS"
dat_hyads$FacID <- str_split(rownames(dat_hyads), ":", simplify=T)[,1] ## get facility IDs
dat_hyads$ZIP <- str_split(rownames(dat_hyads), ":", simplify=T)[,2] ## get zipcodes

dat <- dat_hyads %>% ## merge HyADS with zipcode covariates and facility covariates
  right_join(zip_var, by="ZIP") %>% #ZipcodeData.Rda has more zips?
  full_join(facility_var_2005, by="FacID")
rm("hyads_var")
rm("dat_hyads")

dat$logPop <- log(dat$TotPop)
dat$logMedianHHInc <- log(dat$MedianHHInc)
dat$logPopPerSQM <- log(dat$PopPerSQM)
dat$logHeatInput <- log(dat$totHeatInput)
dat$logOpTime <- log(dat$totOpTime)
#dat$logHyads <- log(dat$tot.hyads)


covariates <- c("logPop", "PctUrban", "PctHisp", "PctHighSchool", "logMedianHHInc", "PctPoor", "PctOccupied", "PctMovedIn5", "logPopPerSQM", #census
                 "smokerate", #smoking
                 "temp", "rh", #weather
                 "totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2", #power plant characteristics
                 "mean_age", "Female_rate", "White_rate", "Black_rate", #Medicare characteristics
                 "HyADS") #HyADS

dat$Y <- dat$IHD_admissions
dat$Z <- ifelse(dat[,'ScrubbedFacility'] == "TRUE", 1, 0)

##### Generate Treatment
X.formula <- paste(covariates, collapse=" + ")
formulaZ <- as.formula(paste("Z ~ ", X.formula))
modZ.log <- glm(formulaZ, family = binomial(), data = dat)
dat$propensity_oracle <- predict(modZ.log, newdata = dat, type = "response")
dat$simZ <- rbinom(nrow(dat),1,dat$propensity_oracle)

##### Generate Outcome
formulaY <- as.formula(paste("Y ~ Z + ", X.formula))
modY.lm <- lm(formulaY, data = dat)
dat <- mutate(dat, simY = -300 -1*Z + 10*logPop + 5*PctUrban - 2*PctHisp - 2*PctHighSchool - 5*logMedianHHInc 
              + 5*PctPoor - 1*PctOccupied - 1*PctMovedIn5 + 1*logPopPerSQM + 10*smokerate + 1*temp + 1*rh 
              + 1*mean_age + 10*Female_rate + 5*White_rate -5*Black_rate + 0.0005*HyADS + rnorm(1,0,10))

##### Propensity Score Analysis #####

## Train-Test Split
train_prop <- 0.75
set.seed(123)

train_ind <- sample(seq_len(nrow(dat)), size = floor(train_prop * nrow(dat)))
df_train <- dat[train_ind, ]
df_test <- dat[-train_ind, ]

train_prop <- df_train['simZ']
train_X <- df_train[-grep('simZ', colnames(df_train))]
test_prop <- df_test['simZ']
test_X <- df_test[-grep('simZ', colnames(df_test))]

xgbMod <- xgboost(data = as.matrix(train_X), label = train_prop, 
                 max_depth = 2,
                 eta = 1,
                 nthread = 2,
                 nrounds = 4,
                 objective = "binary:logistic",
                 verbose = 1)

prop_pred <- predict(xgbMod, test_X)

err <- mean(as.numeric(prop_pred > 0.5) != test_prop)
print(paste("test-error=", err))

