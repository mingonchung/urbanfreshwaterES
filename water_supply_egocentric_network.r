rm(list=ls())

library(lme4)
library(lmerTest)
library(car)
library(memisc)
#library(texreg)
library(MuMIn)

setwd("D:/Min_Gon_Chung/Graduate_PhD/Research/Freshwater_data_030119/7_Data_analysis/process_040921/hlm")

wtr.hlm <- read.csv("freshwater_watersheds_vars_gmis.csv", header=T, stringsAsFactors=FALSE)

### urban population > 300k, 2000-2010
wtr.hlm.sel <- wtr.hlm[which(wtr.hlm$urb_pop_00_10 > 300),]

###########################
### log transformation ####
###########################
### (DV1) water supply
wtr.hlm.sel$lGHM_15_median <- log1p(wtr.hlm.sel$GHM_15_median)

### (DV2) water sediment
wtr.hlm.sel$lsediment_kg_m3 <- log1p(wtr.hlm.sel$sediment_kg_m3)

### (IV1) % of Forest, PAs
wtr.hlm.sel$lforest_pct_PA_basin <- log1p(wtr.hlm.sel$forest_pct_PA_basin)

### (IV1) % of Forest, no PAs
wtr.hlm.sel$lforest_pct_noPA_basin <- log1p(wtr.hlm.sel$forest_pct_noPA_basin)
wtr.hlm.sel$lforest_pct_noPA_basin[is.na(wtr.hlm.sel$lforest_pct_noPA_basin)] <- 0

### (IV1) % of Forest, PAs + no PAs
wtr.hlm.sel$lforest_cover_pct <- log1p(wtr.hlm.sel$forest_cover_pct)

### (IV2) % of Wetland, PAs
wtr.hlm.sel$lwetland_pct_PA_basin <- log1p(wtr.hlm.sel$wetland_pct_PA_basin)

### (IV2) % of Wetland, no PAs
wtr.hlm.sel$lwetland_pct_noPA_basin <- log1p(wtr.hlm.sel$wetland_pct_noPA_basin)
wtr.hlm.sel$lwetland_pct_noPA_basin[is.na(wtr.hlm.sel$lwetland_pct_noPA_basin)] <- 0

### (IV2) % of Wetland, PAs + no PAs
wtr.hlm.sel$lwetland_cover_pct <- log1p(wtr.hlm.sel$wetland_cover_pct)


### (IV3) dam density
wtr.hlm.sel$ldam_density <- log1p(wtr.hlm.sel$dam_density)

### (IV3) dam capacity
wtr.hlm.sel$ldam_capacity <- log1p(wtr.hlm.sel$dam_capacity)

### (IV4) % of irrigation
wtr.hlm.sel$lirrigation_pct <- log1p(wtr.hlm.sel$irrigation_pct)

### (IV5) Watershed areas, km2
wtr.hlm.sel$lbasin_area <- log1p(wtr.hlm.sel$basin_area)

### (IV6) Watershed distance, km
wtr.hlm.sel$ldistance <- log1p(wtr.hlm.sel$distance)

### (IV7) Elevation, meter
#min(wtr.hlm.sel$elevation)
wtr.hlm.sel$elevation_0 <- wtr.hlm.sel$elevation + 3
wtr.hlm.sel$lelevation <- log1p(wtr.hlm.sel$elevation_0)


### (IV8) Slope, degree
wtr.hlm.sel$lslope <- log1p(wtr.hlm.sel$slope)

### (IV9) IWS, 0-1
wtr.hlm.sel$lIWS_program <- log1p(wtr.hlm.sel$IWS_program)

### (IV10) Urban population, 1000 persons
wtr.hlm.sel$lurb_pop_00_10 <- log1p(wtr.hlm.sel$urb_pop_00_10)

### (IV11) Urban GDP per grid, billion USD in 2005 per degree (50km by 50km)
wtr.hlm.sel$lurb_GDP_00_10 <- log1p(wtr.hlm.sel$urb_GDP_00_10)

### (IV12) Temperature, Celsius degree
#min(wtr.hlm.sel$city_temp)
wtr.hlm.sel$city_temp_0 <- wtr.hlm.sel$city_temp + 1.771133
wtr.hlm.sel$lcity_temp <- log1p(wtr.hlm.sel$city_temp_0)

### (IV13) Precipitation, mm
wtr.hlm.sel$lcity_prcp <- log1p(wtr.hlm.sel$city_prcp)

### (IV14) Impervious surface, pct
wtr.hlm.sel$lcity_imper_pct <- log1p(wtr.hlm.sel$city_imper_pct)

##########################
## Model 1: Null model ###
##########################
# (1-1) freshwater
wtr.hlm.mod1 <- lmer(lGHM_15_median ~ 1 + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod1)
rand(wtr.hlm.mod1)


####################################
## Model 2: Watersheds predictors ##
####################################
# (2-1) freshwater, PAs
wtr.hlm.mod2.wtr <- lmer(lGHM_15_median ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod2.wtr)
vif(wtr.hlm.mod2.wtr)
rand(wtr.hlm.mod2.wtr)

# (2-1-1) freshwater, non-PAs
wtr.hlm.mod2.wtr.NPA <- lmer(lGHM_15_median ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod2.wtr.NPA)
vif(wtr.hlm.mod2.wtr.NPA)
rand(wtr.hlm.mod2.wtr.NPA)

# (2-1-2) freshwater, PAs + Non_PAs
wtr.hlm.mod2.wtr.all <- lmer(lGHM_15_median ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod2.wtr.all)
vif(wtr.hlm.mod2.wtr.all)
rand(wtr.hlm.mod2.wtr.all)

# (2-1-3) freshwater, dam capacity instead of dam density
wtr.hlm.mod2.wtr.DC <- lmer(lGHM_15_median ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod2.wtr.DC)
vif(wtr.hlm.mod2.wtr.DC)
rand(wtr.hlm.mod2.wtr.DC)


#############################################
## Model 3: watersheds + urban predictors ###
#############################################
# (3-1) freshwater, PAs
wtr.hlm.mod3.wtr <- lmer(lGHM_15_median ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod3.wtr)
vif(wtr.hlm.mod3.wtr)
rand(wtr.hlm.mod3.wtr)

anova(wtr.hlm.mod2.wtr, wtr.hlm.mod3.wtr)

sink('./results/freshwater_PAs.txt')
print(summary(wtr.hlm.mod3.wtr))

sink('./results/freshwater_PAs_VIF.txt')
print(vif(wtr.hlm.mod3.wtr))

sink()
write.csv(round(summary(wtr.hlm.mod3.wtr)$coefficients,3), './results/freshwater_PAs_coef.csv')

# (3-1-1) freshwater, non-PAs
wtr.hlm.mod3.wtr.NPA <- lmer(lGHM_15_median ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp  + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod3.wtr.NPA)
vif(wtr.hlm.mod3.wtr.NPA)
rand(wtr.hlm.mod3.wtr.NPA)

sink('./results/freshwater_nonPAs.txt')
print(summary(wtr.hlm.mod3.wtr.NPA))
sink()

write.csv(round(summary(wtr.hlm.mod3.wtr.NPA)$coefficients,3), './results/freshwater_nonPAs_coef.csv')


# (3-1-2) freshwater, PAs + Non_PAs
wtr.hlm.mod3.wtr.all <- lmer(lGHM_15_median ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp  + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod3.wtr.all)
vif(wtr.hlm.mod3.wtr.all)
rand(wtr.hlm.mod3.wtr.all)

sink('./results/freshwater_all.txt')
print(summary(wtr.hlm.mod3.wtr.all))
sink()

write.csv(round(summary(wtr.hlm.mod3.wtr.all)$coefficients,3), './results/freshwater_all_coef.csv')

  
  
# (3-1-3) freshwater, dam capacity instead of dam density
wtr.hlm.mod3.wtr.DC <- lmer(lGHM_15_median ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct +  lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp + (1|City_ID), wtr.hlm.sel)
summary(wtr.hlm.mod3.wtr.DC)
vif(wtr.hlm.mod3.wtr.DC)
rand(wtr.hlm.mod3.wtr.DC)

sink('./results/freshwater_PAs_dam_capacity.txt')
print(summary(wtr.hlm.mod3.wtr.DC))
sink()

write.csv(round(summary(wtr.hlm.mod3.wtr.DC)$coefficients,3), './results/freshwater_PAs_dam_capacity_coef.csv')

