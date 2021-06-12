rm(list=ls())

library(lme4)
library(lmerTest)
library(car)
library(memisc)
library(MuMIn)

setwd("D:/Min_Gon_Chung/Graduate_PhD/Research/Freshwater_data_030119/7_Data_analysis/process_040921/hlm")

hydpw.hlm <- read.csv("hydpw_watersheds_vars_gmis.csv", header=T, stringsAsFactors=FALSE)
#hydpw_watersheds_vars.csv

### urban population > 300k, 2000-2010
hydpw.hlm.sel <- hydpw.hlm[which(hydpw.hlm$urb_pop_00_10 > 300),]

###########################
### log transformation ####
###########################
### (DV1) hydropower production
hydpw.hlm.sel$lhydpw_capacity <- log1p(hydpw.hlm.sel$hydpw_capacity)

### (IV1) % of Forest, PAs
hydpw.hlm.sel$lforest_pct_PA_basin <- log1p(hydpw.hlm.sel$forest_pct_PA_basin)

### (IV1) % of Forest, no PAs
hydpw.hlm.sel$lforest_pct_noPA_basin <- log1p(hydpw.hlm.sel$forest_pct_noPA_basin)
hydpw.hlm.sel$lforest_pct_noPA_basin[is.na(hydpw.hlm.sel$lforest_pct_noPA_basin)] <- 0

### (IV1) % of Forest, PAs + no PAs
hydpw.hlm.sel$lforest_cover_pct <- log1p(hydpw.hlm.sel$forest_cover_pct)

### (IV2) % of Wetland, PAs
hydpw.hlm.sel$lwetland_pct_PA_basin <- log1p(hydpw.hlm.sel$wetland_pct_PA_basin)

### (IV2) % of Wetland, no PAs
hydpw.hlm.sel$lwetland_pct_noPA_basin <- log1p(hydpw.hlm.sel$wetland_pct_noPA_basin)
hydpw.hlm.sel$lwetland_pct_noPA_basin[is.na(hydpw.hlm.sel$lwetland_pct_noPA_basin)] <- 0

### (IV2) % of Wetland, PAs + no PAs
hydpw.hlm.sel$lwetland_cover_pct <- log1p(hydpw.hlm.sel$wetland_cover_pct)

### (IV3) dam density
hydpw.hlm.sel$ldam_density <- log1p(hydpw.hlm.sel$dam_density)

### (IV3) dam capacity
hydpw.hlm.sel$ldam_capacity <- log1p(hydpw.hlm.sel$dam_capacity)

### (IV4) % of irrigation
hydpw.hlm.sel$lirrigation_pct <- log1p(hydpw.hlm.sel$irrigation_pct)

### (IV5) Watershed areas, km2
hydpw.hlm.sel$lbasin_area <- log1p(hydpw.hlm.sel$basin_area)

### (IV6) Watershed distance, km
hydpw.hlm.sel$ldistance <- log1p(hydpw.hlm.sel$distance)

### (IV7) Elevation, meter
#min(hydpw.hlm.sel$elevation)
hydpw.hlm.sel$elevation_0 <- hydpw.hlm.sel$elevation
hydpw.hlm.sel$lelevation <- log1p(hydpw.hlm.sel$elevation_0)

### (IV8) Slope, degree
hydpw.hlm.sel$lslope <- log1p(hydpw.hlm.sel$slope)

### (IV9) IWS, 0-1
hydpw.hlm.sel$lIWS_program <- log1p(hydpw.hlm.sel$IWS_program)

### (IV10) Urban population, 1000 persons
hydpw.hlm.sel$lurb_pop_00_10 <- log1p(hydpw.hlm.sel$urb_pop_00_10)

### (IV11) Urban GDP per grid, billion USD in 2005 per degree (50km by 50km)
hydpw.hlm.sel$lurb_GDP_00_10 <- log1p(hydpw.hlm.sel$urb_GDP_00_10)

### (IV12) Temperature, Celsius degree
#min(hydpw.hlm.sel$city_temp)
hydpw.hlm.sel$city_temp_0 <- hydpw.hlm.sel$city_temp
hydpw.hlm.sel$lcity_temp <- log1p(hydpw.hlm.sel$city_temp_0)

### (IV13) Precipitation, mm
hydpw.hlm.sel$lcity_prcp <- log1p(hydpw.hlm.sel$city_prcp)

### (IV14) Impervious surface, pct
hydpw.hlm.sel$lcity_imper_pct <- log1p(hydpw.hlm.sel$city_imper_pct)

##########################
## Model 1: Null model ###
##########################
# (1-1) hydropower
hydpw.hlm.mod1 <- lmer(lhydpw_capacity ~ 1 + (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod1)
rand(hydpw.hlm.mod1)


####################################
## Model 2: Watersheds predictors ##
####################################
# (2-3) hydropower, PAs
hydpw.hlm.mod2 <- lmer(lhydpw_capacity ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod2)
vif(hydpw.hlm.mod2)
rand(hydpw.hlm.mod2)

# (2-3-1) hydropower, non-PAs
hydpw.hlm.mod2.NPA <- lmer(lhydpw_capacity ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod2.NPA)
vif(hydpw.hlm.mod2.NPA)
rand(hydpw.hlm.mod2.NPA)

# (2-3-2) hydropower, PAs + Non_PAs
hydpw.hlm.mod2.all <- lmer(lhydpw_capacity ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod2.all)
vif(hydpw.hlm.mod2.all)
rand(hydpw.hlm.mod2.all)

# (2-3-3) hydropower, dam capacity instead of dam density
hydpw.hlm.mod2.DC <- lmer(lhydpw_capacity ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod2.DC)
vif(hydpw.hlm.mod2.DC)
rand(hydpw.hlm.mod2.DC)


#############################################
## Model 3: watersheds + urban predictors ###
#############################################
# (3-3) hydropower, PAs
hydpw.hlm.mod3 <- lmer(lhydpw_capacity ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod3)
vif(hydpw.hlm.mod3)
rand(hydpw.hlm.mod3)

anova(hydpw.hlm.mod2, hydpw.hlm.mod3)

sink('./results/hydpw_PAs.txt')
print(summary(hydpw.hlm.mod3))

sink('./results/hydpw_PAs_VIF.txt')
print(vif(hydpw.hlm.mod3))

sink()
write.csv(round(summary(hydpw.hlm.mod3)$coefficients,3), './results/hydpw_PAs_coef.csv')


# (3-3-1) hydropower, non-PAs
hydpw.hlm.mod3.NPA <- lmer(lhydpw_capacity ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod3.NPA)
vif(hydpw.hlm.mod3.NPA)
rand(hydpw.hlm.mod3.NPA)

sink('./results/hydpw_nonPAs.txt')
print(summary(hydpw.hlm.mod3.NPA))
sink()

write.csv(round(summary(hydpw.hlm.mod3.NPA)$coefficients,3), './results/hydpw_nonPAs_coef.csv')

# (3-3-2) hydropower, PAs + Non_PAs
hydpw.hlm.mod3.all <- lmer(lhydpw_capacity ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp  + (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod3.all)
vif(hydpw.hlm.mod3.all)
rand(hydpw.hlm.mod3.all)

sink('./results/hydpw_all.txt')
print(summary(hydpw.hlm.mod3.all))
sink()

write.csv(round(summary(hydpw.hlm.mod3.all)$coefficients,3), './results/hydpw_all_coef.csv')

# (3-3-3) hydropower, dam capacity instead of dam density
hydpw.hlm.mod3.DC <- lmer(lhydpw_capacity ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), hydpw.hlm.sel)
summary(hydpw.hlm.mod3.DC)
vif(hydpw.hlm.mod3.DC)
rand(hydpw.hlm.mod3.DC)

sink('./results/hydpw_PAs_dam_capacity.txt')
print(summary(hydpw.hlm.mod3.DC))
sink()

write.csv(round(summary(hydpw.hlm.mod3.DC)$coefficients,3), './results/hydpw_PAs_coef_dam_capacity.csv')
