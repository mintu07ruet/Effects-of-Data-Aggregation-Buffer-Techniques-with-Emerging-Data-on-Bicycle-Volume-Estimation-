rm(list = ls())

if (!require("pacman")) install.packages("pacman")
p_load(dplyr)
#p_load(MASS)
p_load(sandwich)
#library(mgcv)
p_load(lmtest)
p_load(caret)
p_load(AER)
p_load(cvTools)
p_load(ggplot2)
p_load(readr)
#library(plyr)
#library(glmmTMB)

##### FUNCTION DEFINITIONS #####
cvstats <- function (p, k, r=1) {
  n = k * r
  stats <- p %>%
    group_by(Resample) %>%
    summarize(rmse=sqrt(mean((pred - obs) ^ 2)),
              mae=mean(abs(pred - obs)),
              mape=mean(abs((pred - obs) / obs * 100)),
              r2=cor(pred, obs)^2,
              mfr2=((sum()))) %>%
    summarize(mean(rmse), se_rmse=sd(rmse)/sqrt(n), mean(mae), 
              se_mae=sd(mae)/sqrt(n), mean(mape), se_mape=sd(mape)/sqrt(n),
              mean(r2), se_r2=sd(r2)/sqrt(n))
  stats
}

cvstats.tertile <- function(p, k, r=1) {
  # Note: assignment is stochastic at boundaries, thus the fractional n's
  n_cv = k * r  # total out of sample predictions
  
  stats <- p %>% 
    mutate(tertile = ntile(obs, 3)) %>%
    group_by(Resample, tertile) %>%
    summarize(rmse=sqrt(mean((pred - obs) ^ 2)),
              mae=mean(abs(pred - obs)),
              mape=mean(abs((pred - obs) / obs * 100)),
              r2=cor(pred, obs, use="complete.obs")^2, n=n()) %>%
    group_by(tertile) %>%
    summarize(mean(rmse), se_rmse=sd(rmse)/sqrt(n_cv), mean(mae), 
              se_mae=sd(mae)/sqrt(n_cv), mean(mape), se_mape=sd(mape)/sqrt(n_cv),
              mean(r2), se_r2=sd(r2)/sqrt(n_cv), n=sum(n)/r)
  
  stats
}

cvstats.region <- function(m, d, k, r) {
  # for now group must be in select var list
  #n_cv = k * r
  
  p2 <- d %>% select(rowIndex, site_id, region, AADBT, site_name) %>%
    inner_join(m$pred, by="rowIndex") 
  
  n_reg <- p2 %>%
    distinct(rowIndex, region) %>%
    group_by(region) %>%
    summarize(n=n())
  
  stats <- p2 %>% 
    group_by(Resample, region) %>%
    summarize(rmse=sqrt(mean((pred - obs) ^ 2)),
              mae=mean(abs(pred - obs)),
              mape=mean(abs((pred - obs) / obs * 100)),
              ae=mean(pred-obs),
              ape=mean((pred - obs) / obs * 100),
              r2=cor(pred, obs, use="complete.obs")^2,
              n=n()) %>%
    group_by(region) %>%
    summarize(mean(rmse), mean(mae), 
              mean(mape), mean(r2, na.rm=T)) %>%
    left_join(n_reg, by="region")
  
  stats
}

cvstats.pattern <- function(m, d, k, r) {
  # for now group must be in select var list
  p2 <- d %>% select(rowIndex, site_id, AADBT, site_name, travel_pattern) %>%
    inner_join(m$pred, by="rowIndex")
  
  stats <- p2 %>% 
    group_by(travel_pattern) %>%
    summarize(rmse=sqrt(mean((pred - obs) ^ 2)),
              mae=mean(abs(pred - obs)),
              mape=mean(abs((pred - obs) / obs * 100)),
              avgpe=mean((pred - obs) / obs * 100),
              avgerr=mean(pred-obs))
  print(stats)
}

mf_r2 <- function(m) {
  # McFadden's pseudo-R2
  r2 <- 1 - (m$deviance / m$null.deviance)
  r2
}

summarize_model <- function(m, d, k, r, robust_se="HC3") {
  # m - model train/test results from caret package
  # k - number of Cv folds
  # r - number of fold repeats
  p <- m$pred
  print(cvstats(p, k, r))
  print(cvstats.tertile(p, k, r))
  print(cvstats.region(m, d, k, r))
  print(summary(m))
  print("")
  print(paste("McFadden's pseudo-R2 =", mf_r2(m$finalModel)))
  print("")
  print(coeftest(m$finalModel, 
                 vcov = vcovHC(m$finalModel, type=robust_se)))  # robust SEs 
}

plot_regions <- function(m, d, title="", scales="fixed") {
  p2 <- d %>% select(rowIndex, site_id, region, AADBT, site_name) %>%
    inner_join(m$pred, by="rowIndex") %>%
    group_by(site_id, region, site_name) %>%
    summarize(avg_obs=mean(obs), avg_pred=mean(pred))
  
  p2 %>% ggplot(aes(x=avg_pred, y=avg_obs)) +
    geom_point(col="blue") +
    geom_abline(intercept=0, slope=1, lty=2) +
    facet_wrap(~ region, scales=scales) +
    labs(x="mean predicted AADB", y="observed AADB") +
    ggtitle(title)
}

plot_all <- function(m, d, title="", label_outliers=0) {
  # label_outliers - label if avg err falls in top X pct of cases
  p2 <- d %>% select(rowIndex, site_id, region, AADBT, site_name) %>%
    inner_join(m$pred, by="rowIndex") %>%
    group_by(site_id, region, site_name) %>%
    summarize(avg_obs=mean(obs), avg_pred=mean(pred)) %>%
    mutate(avg_err=mean(avg_pred - avg_obs)) 
  label_min <- quantile(abs(p2$avg_err), p=((100 - label_outliers) / 100))
  
  p2 %>% ggplot(aes(x=avg_pred, y=avg_obs, col=as.factor(region),
                    shape=as.factor(region))) +
    geom_abline(intercept=0, slope=1, lty=2) +
    geom_point(size=2.0, fill=NA) +
    labs(x="mean predicted AADB", y="observed AADB") + 
    geom_text(data=subset(p2, avg_err >= 0 & abs(avg_err) >= label_min),
              aes(avg_pred, avg_obs, label=substr(site_name, 1, 11)),
              nudge_y=30, size=3) +
    geom_text(data=subset(p2, avg_err < 0 & abs(avg_err) >= label_min),
              aes(avg_pred, avg_obs, label=substr(site_name, 1, 11)),
              nudge_y=-30, size=3) +
    ggtitle(title)
}

##### END FUNCTION DEFINITIONS #####

# Import Data
all_data <- read.csv("E:/Bike Fusion/New Plan/Circular Automatic Collected data for modeling/Buffer Paper/R_exported_buffer_paper_2019_Six_cities_data.csv",stringsAsFactors = F)

# Get the OSM facilities at each count station by counter types
facilities=subset(all_data, select=c("type","region","primary_binary","secondary_binary","tertiary_binary","residential_binary","path_binary","cycleway_binary","cycleway_lane_binary","cycleway_track_all_binary","footway_binary"))
#require(dplyr)
# count permanent and short term sites
facilities %>% count(facilities$type,facilities$region)
# sum  facilities by type and region

lapply(facilities, count)
pc=subset(facilities,type=='permanent')
lapply(pc, count)
# count facilities
pc1=subset(pc, select=c("region","primary_binary","secondary_binary","tertiary_binary","residential_binary","path_binary","cycleway_binary","cycleway_lane_binary","cycleway_track_all_binary","footway_binary"))
require(data.table)
dt1 <- data.table(pc1)
dt1.sum <- dt1[, lapply(.SD, sum), by = 'region']
dt1.sum




sc=subset(facilities,type=='short-term')
lapply(sc, count)
# count facilities
sc1=subset(sc, select=c("region","primary_binary","secondary_binary","tertiary_binary","residential_binary","path_binary","cycleway_binary","cycleway_lane_binary","cycleway_track_all_binary","footway_binary"))
require(data.table)
dt <- data.table(sc1)
dt.sum <- dt[, lapply(.SD, sum), by = 'region']
dt.sum




# Run Air+Network Model for all buffers for six cities
# Set up train/test controller 
seed = 37914 #37914
k = 10  # folds
r = 5  # repeats
set.seed(seed)  # note checked fairly stable w/ seed changes

# Note: eventually went back to caret w/ stratified folds 
train_control5 <- trainControl(method="repeatedcv", number=k, repeats=r,
                               savePredictions = TRUE)



# Run Euclidian buffer model
# Buffer Size =all
# run the model
quantile(as.integer(all_data$AADBT), p=seq(0, 1, 1/3))

set.seed(seed)
pm8.cv <- train(as.integer(AADBT) ~   stv_c_adb
                +                log_stl_raw
                +           Tertiary_ohm
                +      Student.Access_tm
                +     Number.of.jobs_tm
                +        log_stv_nc_adb
                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm8.cv, all_data, k, r)
plot_regions(pm8.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm8.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)  # label outliers in top X%, NA turns off labeling


# 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~   stv_c_adb
                +                        tertiary_binary
                +                     Distance.to.forest
                +                             stv_nc_adb
                +                            log_stv_stl
                +                           Footway_otm
                +    pct_at_least_college_education_otm
                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, all_data, k, r)
plot_regions(pm1.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)  # label outliers in top X%, NA turns off labeling

# 0.25 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~  stv_c_adb
                +                       tertiary_binary
                +                    Distance_to_CBD_mi
                +                            stv_nc_adb
                +                           log_stv_stl
                +                         Bus.Stops_qm
                +                           Footway_qm

            
,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, all_data, k, r)
plot_regions(pm2.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)


# 0.5 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~      stv_c_adb
                +                       tertiary_binary
                +                    Distance.to.forest
                +                            stv_nc_adb
                +                           log_stv_stl
                +                    Student.Access_hm
                +    pct_at_least_college_education_hm
                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, all_data, k, r)
plot_regions(pm3.cv, all_data, scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, all_data, title="Pool Model for 0.5 mile Euclidian and Network Buffers", label_outliers=NA)


# 0.75 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~           stv_c_adb
                +                                log_stl_raw
                +                        tertiary_binary
                +                     Distance.to.forest
                +                             stv_nc_adb
                +                     Distance_to_CBD_mi
                +    pct_at_least_college_education_tfm

                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, all_data, k, r)
plot_regions(pm4.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)

# 1.0 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~             stv_c_adb
                +                   min_dist_to_college
                +                            stv_nc_adb
                +                           log_stv_stl
                +                    Student.Access_om

                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, all_data, k, r)
plot_regions(pm5.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)

# 1.5 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~             stv_c_adb
                +                     Student.Access_ohm
                +                           log_stv_stl
                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, all_data, k, r)
plot_regions(pm6.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)

# 2.0 mile buffer model
set.seed(seed)
pm7.cv <- train(as.integer(AADBT) ~             stv_c_adb
                +            Student.Access_tm
                +                   stv_nc_adb
                +                 log_stv_stl
                ,
                data = all_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm7.cv, all_data, k, r)
plot_regions(pm7.cv, all_data, title="Pool Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm7.cv, all_data, title="Pool Model for all Buffers", label_outliers=NA)


##  Portland  Model
p_data <- all_data %>% filter(region %in% c("Portland"))
p_data$rowIndex <- as.integer(rownames(p_data))
quantile(as.integer(p_data$AADBT), p=seq(0, 1, 1/3))
# all mile buffer model
set.seed(seed)
pm0.cv <- train(as.integer(AADBT) ~     log_stv_c_adb
                +                University_tm
                +Cycleway_tm
                +Distance_to_CBD_mi
              			
                               ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm0.cv, p_data, k, r)
plot_regions(pm0.cv, p_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm0.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 


# Portland 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~      tertiary_binary
                +                      path_binary
                +     Distance.to.Residential.Area
                +                    log_stv_c_adb
                +              Distance_to_CBD_mi
                +             Industrial.Area_otm
                +            Residential_Road_otm
                +                  Water.Area_otm
                +                   Secondary_otm
                +                  HH_density_otm
                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, p_data, k, r)
plot_regions(pm1.cv, p_data, title="Portland Best Model with 0.10-mile Air Buffer", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 


# Portland 0.25 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~         tertiary_binary
                +     Distance.to.Residential.Center
                +                      log_stv_c_adb
                +                           Path_qm

                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, p_data, k, r)
plot_regions(pm2.cv, p_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Portland 0.5 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~         stv_c_adb
                +                       tertiary_binary
                +                   Distance_to_CBD_mi
                +Distance.to.Industrial.Area
                +pct_at_least_college_education_hm
	
                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, p_data, k, r)
plot_regions(pm3.cv, p_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Portland 0.75 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~         tertiary_binary
                +        Distance.to.Industrial.Area
                +     Distance.to.Residential.Center
                +                      log_stv_c_adb
                +                 Distance_to_CBD_mi
                +              Residential_Road_tfm



                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, p_data, k, r)
plot_regions(pm4.cv, p_data, title="Portland Best Model with 0.75 Mile Air and Network Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, p_data, title="Portland Best Model with 0.75 Mile Buffers", label_outliers=NA) 

# Portland 1.0 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~   stv_c_adb
                +                       tertiary_binary
#                +                           path_binary
#                +         Distance.to.Industrial.Center
                +        Distance.to.Residential.Center
#                +              Distance.to.Grass.Center
#                +                          Point.Bridge
#                +                        BikeFac_binary
#                +               Intersection_Density_om
                +                          Cycleway_om
#                +                cycleway_track_all_om
                +                          pct_male_om
#                +    pct_at_least_college_education_om
#                +               Industrial.Area_om_net
#                +                   Retail.Area_om_net
#                +                    Grass.Area_om_net
#                +                     Secondary_om_net
#                +                Student.Access_om_net
                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, p_data, k, r)
plot_regions(pm5.cv, p_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Portland 1.5 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~   stv_c_adb
                +                         tertiary_binary
#                +                             path_binary
#                +      Distance.to.Commercial.Area.Center
                +            Distance.to.Residential.Area
#                +                 Distance.to.Retail.Area
#                +                Distance.to.Grass.Center
#                +                     min_dist_to_college
#                +                             Point.Speed
#                +                    Industrial.Area_ohm
#                +                        Retail.Area_ohm
#                +                           Cycleway_ohm
#                +                           pct_male_ohm
#                +    Percentage.of.Bike.Commuter_ohm_net
#                +               Industrial.Area_ohm_net
#                +                    Grass.Area_ohm_net
#                +                        Footway_ohm_net
#                +                 Student.Access_ohm_net
#                +                         slope_ohm_net
#                +                     sep_bikeway_binary
#                +                        sep_bikeway_ohm
                +                     Distance_to_CBD_mi
                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, p_data, k, r)
plot_regions(pm6.cv, p_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Portland 2.0 mile buffer model
set.seed(seed)
pm7.cv <- train(as.integer(AADBT) ~    stv_c_adb
                +                        tertiary_binary
                +                     Industrial.Area_tm
                +                         University_tm
                +                           Cycleway_tm

                ,
                data = p_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm7.cv, p_data, k, r)
plot_regions(pm7.cv, p_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm7.cv, p_data, title="Portland Model for all Buffers", label_outliers=NA) 


# Lets Run with Eugene Model
e_data <- all_data %>% filter(region %in% c("Eugene"))
e_data$rowIndex <- as.integer(rownames(e_data))
quantile(as.integer(e_data$AADBT), p=seq(0, 1, 1/3))
# all mile buffer model
set.seed(seed)
pm0.cv <- train(as.integer(AADBT) ~ Bus.Stops_qm
                +                          Footway_hm
                +                       pct_white_tfm
                +                          college_tm
                +                Median_HH_income_tm
                +                      log_stv_c_adb
                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm0.cv, e_data, k, r)
plot_regions(pm0.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm0.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Eugene 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~  Distance.to.Retail.Center
                +                       Distance.to.Park
                +                             stv_nc_adb
                +                       Forest.Area_otm
                +                   Bicycle.Parking_otm
                +                       University_otm
                +                  Median_HH_income_otm

                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, e_data, k, r)
plot_regions(pm1.cv, e_data, title="Eugene Best Model with 0.1 mile Air and Network Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, e_data, title="Eugene Best Model with 0.1 mile Air and Network Buffers", label_outliers=NA) 

# Eugene 0.25 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~  stl_raw
                +                        log_stv_nc_adb
                +                        University_qm
                +                    Student.Access_qm
                +    pct_at_least_college_education_qm

                
                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, e_data, k, r)
plot_regions(pm2.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Eugene 0.50 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~ Distance.to.Retail.Center
                +                   Distance.to.Park.Center
                +                            log_stv_nc_adb
                +                              log_stl_raw
                +                              Bus.Stops_hm
               
                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, e_data, k, r)
plot_regions(pm3.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Eugene 0.75 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~   Distance.to.Park
                +               Distance.to.Water.Center
                +                         log_stv_nc_adb
                +                            log_stl_raw
                +                               Path_tfm
                +                population_density_tfm
                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, e_data, k, r)
plot_regions(pm4.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Eugene 1.0 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~   log_stv_c_adb
                +                       tertiary_binary
                +                Distance.to.Water.Body
                +                        log_stv_nc_adb
                +                          maxspeed_om
                +    pct_at_least_college_education_om


                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, e_data, k, r)
plot_regions(pm4.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Eugene 1.5 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~   Distance.to.Retail.Area
                +                           Distance.to.Park
                +                        Bicycle.Parking_ohm
                +                            log_stv_nc_adb
                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, e_data, k, r)
plot_regions(pm5.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 

# Eugene 2.0 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~  Distance.to.Park.Center
                +      Distance.to.Water.Center
                +                 University_tm
                +               log_stv_nc_adb
                +                 log_stl_raw
                ,
                data = e_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, e_data, k, r)
plot_regions(pm6.cv, e_data, title="Portland Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, e_data, title="Portland Model for all Buffers", label_outliers=NA) 



# Bend Model 
# All buffers
be_data <- all_data %>% filter(region %in% c("Bend"))
be_data$rowIndex <- as.integer(rownames(be_data))
quantile(as.integer(be_data$AADBT), p=seq(0, 1, 1/3))
# Bend all mile buffer model
set.seed(seed)
pm0.cv <- train(as.integer(AADBT) ~  Bus.Stops_qm
                +                  sep_onstreet_binary
                +                        log_stv_c_adb

                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm0.cv, be_data, k, r)
plot_regions(pm0.cv, be_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm0.cv, be_data, title="Bend Model for all Buffers", label_outliers=NA)


# Bend 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                   sep_onstreet_binary
                +                           log_stl_raw
          
                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, be_data, k, r)
plot_regions(pm1.cv, be_data, title="Bend Model for 0.1 mile  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, be_data, title="Bend Model for 0.1 mile Buffers", label_outliers=NA)

# Bend 0.25 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                    sep_onstreet_binary
                +                           log_stl_raw
                +                          Bus.Stops_qm
                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, be_data, k, r)
plot_regions(pm2.cv, be_data, title="Bend Best Model with 0.25-mile Air Buffer", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, be_data, title="Bend Model for 0.25 mile Buffers", label_outliers=NA)

# Bend 0.50 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                   sep_onstreet_binary
                +                           log_stl_raw
                +                    Bicycle.Parking_hm

                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, be_data, k, r)
plot_regions(pm3.cv, be_data, title="Bend Model for 0.25 mile  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, be_data, title="Bend Model for 0.25 mile Buffers", label_outliers=NA)

# Bend 0.750 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~cycleway_binary
                +                       Distance.to.Park
                +                          log_stv_c_adb
                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, be_data, k, r)
plot_regions(pm4.cv, be_data, title="Bend Model for 0.75 mile  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, be_data, title="Bend Model for 0.75 mile Buffers", label_outliers=NA)

# Bend 1.0 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~pct_African_American_om
                +           Distance.to.Grass.Center
                +                sep_onstreet_binary
                +                     log_stv_c_adb
                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, be_data, k, r)
plot_regions(pm5.cv, be_data, title="Bend Model for 0.75 mile  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, be_data, title="Bend Model for 0.75 mile Buffers", label_outliers=NA)

# Bend 1.5 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                        sep_onstreet_binary
                +                               log_stl_raw
                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, be_data, k, r)
plot_regions(pm6.cv, be_data, title="Bend Model for 1.5 mile  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, be_data, title="Bend Model for 1.5 mile Buffers", label_outliers=NA)

# Bend 2.0 mile buffer model
set.seed(seed)
pm7.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +           sep_onstreet_binary
                +                   log_stl_raw
                ,
                data = be_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm7.cv, be_data, k, r)
plot_regions(pm7.cv, be_data, title="Bend Model for 2. mile  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm7.cv, be_data, title="Bend Model for 2.0 mile Buffers", label_outliers=NA)




# Dallas all buffer model
da_data <- all_data %>% filter(region %in% c("Dallas"))
da_data$rowIndex <- as.integer(rownames(da_data))
quantile(as.integer(da_data$AADBT), p=seq(0, 1, 1/3))
# Dallas all mile buffer model
set.seed(seed)
pm0.cv <- train(as.integer(AADBT) ~pct_white_tfm
                +    pct_at_least_college_education_tm
                +                           stv_nc_adb
                +                          log_stl_raw
                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm0.cv, da_data, k, r)
plot_regions(pm0.cv, da_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm0.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)



# Dallas 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~pct_white_otm+Distance_to_Water_Body_mi
                +                                stv_nc_adb
                +                            BikeFac_binary
                +                               log_stl_raw
               ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, da_data, k, r)
plot_regions(pm1.cv, da_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)


# Dallas 0.25 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                            log_stl_raw
                +                    Commercial.Area_qm
                +                             School_qm
                +     pct_at_least_college_education_qm

                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, da_data, k, r)
plot_regions(pm2.cv, da_data, title="Dallas Best Model with 0.25-mile Air Buffer", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)

# Dallas 0.50 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                           log_stl_raw
                
                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, da_data, k, r)
plot_regions(pm3.cv, da_data, title="Dallas Best Model with 0.5 mile Air and Network Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, da_data, title="Dallas Best Model with 0.5 mile Air and Network Buffers", label_outliers=NA)

# Dallas 0.75 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                            log_stl_raw+BikeFac_binary
#                +                        Forest.Area_tfm
#                +               pct_African_American_tfm
#                +              Residential_Road_tfm_net
                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, da_data, k, r)
plot_regions(pm4.cv, da_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)

# Dallas 1.0 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~stv_nc_adb
                +                        BikeFac_binary
                +                           log_stl_raw
                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, da_data, k, r)
plot_regions(pm5.cv, da_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)

# Dallas 1.5 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~
#                  pct_African_American_ohm
#                +                Median.Age_ohm
#                +       Industrial.Area_ohm_net
#                +              Tertiary_ohm_net
                                   stv_nc_adb
                +               BikeFac_binary
                +                  log_stl_raw
                +    Distance_to_Water_Body_mi
                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, da_data, k, r)
plot_regions(pm6.cv, da_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)

# Dallas 2.0 mile buffer model
set.seed(seed)
pm7.cv <- train(as.integer(AADBT) ~pct_female_tm
#                +                             Median.Age_tm
 #               +                    Industrial.Area_tm_net
                +                               stv_nc_adb
                +                              log_stl_raw
                ,
                data = da_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm7.cv, da_data, k, r)
plot_regions(pm7.cv, da_data, title="Bend Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm7.cv, da_data, title="Bend Model for all Buffers", label_outliers=NA)



# Boulder all buffer model
bo_data <- all_data %>% filter(region %in% c("Boulder"))
bo_data$rowIndex <- as.integer(rownames(bo_data))
quantile(as.integer(bo_data$AADBT), p=seq(0, 1, 1/3))
# Boulder all mile buffer model
set.seed(seed)
pm0.cv <- train(as.integer(AADBT) ~pct_at_least_college_education_tm
                +                        Median.Age_tm
                +                        log_stv_c_adb

                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm0.cv, bo_data, k, r)
plot_regions(pm0.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm0.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~min_dist_to_college
                +              Point.Speed
                +             pct_male_otm
                +          log_stv_nc_adb
                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, bo_data, k, r)
plot_regions(pm1.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 0.25 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~log_stv_nc_adb
               +             Median_HH_income_qm
                +               Number.of.jobs_qm
                +                        slope_qm
        #        +     Intersection_Density_qm_net
        #        +         Residential_Road_qm_net
        #        +              Retail.Area_qm_net
       #         +                   School_qm_net
       #         +                 Cycleway_qm_net
      #          +                  bridge_qm_net
      #          +    pct_African_American_qm_net
      #          +          Student.Access_qm_net
       #        +              Median.Age_qm_net
        #        +                     BikeFac_qm
       #         +                  Park_acres_qm
                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, bo_data, k, r)
plot_regions(pm1.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 0.50 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~log_stv_nc_adb
                +                   Residential_Road_hm
                +                        Forest.Area_hm
                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, bo_data, k, r)
plot_regions(pm2.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 0.750 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~log_stv_nc_adb
#                +                         Water.Area_tfm
#                +                          Secondary_tfm
#                +    pct_at_least_college_education_tfm
                +                  Median_HH_income_tfm
#                +                             slope_tfm
                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, bo_data, k, r)
plot_regions(pm3.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 1.0 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~ stv_nc_adb+
#                +                         stv_c_adb
                +       Distance.to.Commercial.Area
                +    Distance.to.Industrial.Center
                +          Distance.to.Retail.Area
#                +               Distance.to.forest
                +                       Number.of.jobs_om_net
                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, bo_data, k, r)
plot_regions(pm4.cv, bo_data, title="Boulder Best Model with 1.0-mile Air Buffer", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 1.5 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~Median_HH_income_ohm
                +                        log_stv_nc_adb
                +                        Park_acres_ohm

                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, bo_data, k, r)
plot_regions(pm5.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Boulder 2.0 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~log_stv_nc_adb
                +Tertiary_tm
                + maxspeed_tm
#                +Distance_to_CBD_mi
                
                ,
                data = bo_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, bo_data, k, r)
plot_regions(pm6.cv, bo_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, bo_data, title="Boulder Model for all Buffers", label_outliers=NA)


# Charlotte all buffer model
ch_data <- all_data %>% filter(region %in% c("Charlotte"))
ch_data$rowIndex <- as.integer(rownames(ch_data))
quantile(as.integer(ch_data$AADBT), p=seq(0, 1, 1/3))
# Charlotte all mile buffer model
set.seed(seed)
pm0.cv <- train(as.integer(AADBT) ~BikeFac_binary
                +              arterial_binary

                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm0.cv, ch_data, k, r)
plot_regions(pm0.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm0.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Charlotte 0.1 mile buffer model
set.seed(seed)
pm1.cv <- train(as.integer(AADBT) ~arterial_binary
                +               Park_acres_otm

                
                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm1.cv, ch_data, k, r)
plot_regions(pm1.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm1.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Charlotte 0.25 mile buffer model
set.seed(seed)
pm2.cv <- train(as.integer(AADBT) ~Median.Age_qm
                +              Park_acres_qm
                
                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm2.cv, ch_data, k, r)
plot_regions(pm2.cv, ch_data, title="Charlotte Best Model with 0.25-mile Air Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm2.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

# Charlotte 0.50 mile buffer model
set.seed(seed)
pm3.cv <- train(as.integer(AADBT) ~Park_acres_hm
                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm3.cv, ch_data, k, r)
plot_regions(pm3.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm3.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

# # Charlotte 0.75 mile buffer model
set.seed(seed)
pm4.cv <- train(as.integer(AADBT) ~Park_acres_tfm
                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm4.cv, ch_data, k, r)
plot_regions(pm4.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm4.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

# # Charlotte 1.0 mile buffer model
set.seed(seed)
pm5.cv <- train(as.integer(AADBT) ~BikeFac_binary
                +              arterial_binary
                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm5.cv, ch_data, k, r)
plot_regions(pm5.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm5.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)


# # Charlotte 1.5 mile buffer model
set.seed(seed)
pm6.cv <- train(as.integer(AADBT) ~BikeFac_binary
                +              arterial_binary
                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm6.cv, ch_data, k, r)
plot_regions(pm6.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm6.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

# # Charlotte 1.5 mile buffer model
set.seed(seed)
pm7.cv <- train(as.integer(AADBT) ~BikeFac_binary
                +              arterial_binary

                ,
                data = ch_data, method = "glm",family='poisson',
                trControl = train_control5)
#family=sm.families.NegativeBinomial(link=sm.families.links.identity)
#family='poisson'
# for nargative Binomial: family="nbinom2"
summarize_model(pm7.cv, ch_data, k, r)
plot_regions(pm7.cv, ch_data, title="Boulder Model for all  Buffers", scales="fixed")  # scales="free" is an option
plot_all(pm7.cv, ch_data, title="Boulder Model for all Buffers", label_outliers=NA)

