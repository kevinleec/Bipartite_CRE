##### Load Libraries & Data #####
library(gdata)
library(stringr)
library(dplyr)
library(data.table)
library(xgboost)

load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/HyADSmat.Rda")
load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/out.zip_pp.rda")
load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/facilities_for_analysis.Rda")
load("~/shared_space/ci3_analysis/zigler_lab/projects/BipartiteInterference_GPS/BipartiteInterference_GPS/Data/ZipcodeData.Rda")

hyads_var <- sumMAP.2005.fac
zip_var <- data
facility_var_2005 <- dat_facility_2005
rm("sumMAP.2005.fac")


##### Convert data to proper format #####
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

X.formula <- paste(covariates, collapse=" + ")
formulaZ <- as.formula(paste("Z ~ ", X.formula))
modZ.log <- glm(formulaZ, family = binomial(), data = dat)
dat$propensity_oracle <- predict(modZ.log, newdata = dat, type = "response")
dat$simZ <- rbinom(nrow(dat),1,dat$propensity_oracle)

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
                 max_depth = 3,
                 eta = 1,
                 nthread = 2,
                 nrounds = 4,
                 objective = "binary:logistic",
                 verbose = 1)

prop_pred <- predict(xgbMod, test_X)

err <- mean(as.numeric(prop_pred > 0.5) != test_prop)
print(paste("test-error=", err))

