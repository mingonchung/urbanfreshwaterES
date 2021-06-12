rm(list=ls())

library(lme4)
library(lmerTest)
library(car)
library(memisc)
#library(texreg)
library(MuMIn)

setwd("D:/Min_Gon_Chung/Graduate_PhD/Research/Freshwater_data_030119/7_Data_analysis/process_040921/hlm")

sed.hlm <- read.csv("freshwater_watersheds_vars_gmis.csv", header=T, stringsAsFactors=FALSE)

### urban population > 300k, 2000-2010
sed.hlm.sel <- sed.hlm[which(sed.hlm$urb_pop_00_10 > 300),]

###########################
### log transformation ####
###########################
### (DV1) water supply
sed.hlm.sel$lGHM_15_median <- log1p(sed.hlm.sel$GHM_15_median)

### (DV2) water sediment
sed.hlm.sel$lsediment_kg_m3 <- log1p(sed.hlm.sel$sediment_kg_m3)

### (IV1) % of Forest, PAs
sed.hlm.sel$lforest_pct_PA_basin <- log1p(sed.hlm.sel$forest_pct_PA_basin)

### (IV1) % of Forest, no PAs
sed.hlm.sel$lforest_pct_noPA_basin <- log1p(sed.hlm.sel$forest_pct_noPA_basin)
sed.hlm.sel$lforest_pct_noPA_basin[is.na(sed.hlm.sel$lforest_pct_noPA_basin)] <- 0

### (IV1) % of Forest, PAs + no PAs
sed.hlm.sel$lforest_cover_pct <- log1p(sed.hlm.sel$forest_cover_pct)

### (IV2) % of Wetland, PAs
sed.hlm.sel$lwetland_pct_PA_basin <- log1p(sed.hlm.sel$wetland_pct_PA_basin)

### (IV2) % of Wetland, no PAs
sed.hlm.sel$lwetland_pct_noPA_basin <- log1p(sed.hlm.sel$wetland_pct_noPA_basin)
sed.hlm.sel$lwetland_pct_noPA_basin[is.na(sed.hlm.sel$lwetland_pct_noPA_basin)] <- 0

### (IV2) % of Wetland, PAs + no PAs
sed.hlm.sel$lwetland_cover_pct <- log1p(sed.hlm.sel$wetland_cover_pct)


### (IV3) dam density
sed.hlm.sel$ldam_density <- log1p(sed.hlm.sel$dam_density)

### (IV3) dam capacity
sed.hlm.sel$ldam_capacity <- log1p(sed.hlm.sel$dam_capacity)

### (IV4) % of irrigation
sed.hlm.sel$lirrigation_pct <- log1p(sed.hlm.sel$irrigation_pct)

### (IV5) Watershed areas, km2
sed.hlm.sel$lbasin_area <- log1p(sed.hlm.sel$basin_area)

### (IV6) Watershed distance, km
sed.hlm.sel$ldistance <- log1p(sed.hlm.sel$distance)

### (IV7) Elevation, meter
#min(sed.hlm.sel$elevation)
sed.hlm.sel$elevation_0 <- sed.hlm.sel$elevation + 3
sed.hlm.sel$lelevation <- log1p(sed.hlm.sel$elevation_0)


### (IV8) Slope, degree
sed.hlm.sel$lslope <- log1p(sed.hlm.sel$slope)

### (IV9) IWS, 0-1
sed.hlm.sel$lIWS_program <- log1p(sed.hlm.sel$IWS_program)

### (IV10) Urban population, 1000 persons
sed.hlm.sel$lurb_pop_00_10 <- log1p(sed.hlm.sel$urb_pop_00_10)

### (IV11) Urban GDP per grid, billion USD in 2005 per degree (50km by 50km)
sed.hlm.sel$lurb_GDP_00_10 <- log1p(sed.hlm.sel$urb_GDP_00_10)

### (IV12) Temperature, Celsius degree
#min(sed.hlm.sel$city_temp)
sed.hlm.sel$city_temp_0 <- sed.hlm.sel$city_temp + 1.771133
sed.hlm.sel$lcity_temp <- log1p(sed.hlm.sel$city_temp_0)

### (IV13) Precipitation, mm
sed.hlm.sel$lcity_prcp <- log1p(sed.hlm.sel$city_prcp)

### (IV14) Impervious surface, pct
sed.hlm.sel$lcity_imper_pct <- log1p(sed.hlm.sel$city_imper_pct)

##########################
## Model 1: Null model ###
##########################
# (1-2) sediment
sed.hlm.mod1 <- lmer(lsediment_kg_m3 ~ 1 + (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod1)
rand(sed.hlm.mod1)

####################################
## Model 2: Watersheds predictors ##
####################################

# (2-2) sediment, PAs
sed.hlm.mod2.sed <- lmer(lsediment_kg_m3 ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod2.sed)
vif(sed.hlm.mod2.sed)
rand(sed.hlm.mod2.sed)

# (2-2-1) sediment, non-PAs
sed.hlm.mod2.sed.NPA <- lmer(lsediment_kg_m3 ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod2.sed.NPA)
vif(sed.hlm.mod2.sed.NPA)
rand(sed.hlm.mod2.sed.NPA)

# (2-2-2) sediment, PAs + Non_PAs
sed.hlm.mod2.sed.all <- lmer(lsediment_kg_m3 ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod2.sed.all)
vif(sed.hlm.mod2.sed.all)
rand(sed.hlm.mod2.sed.all)

# (2-2-3) sediment, dam capacity instead of dam density
sed.hlm.mod2.sed.DC <- lmer(lsediment_kg_m3 ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod2.sed.DC)
vif(sed.hlm.mod2.sed.DC)
rand(sed.hlm.mod2.sed.DC)


#############################################
## Model 3: watersheds + urban predictors ###
#############################################
# (3-2) sediment, PAs
sed.hlm.mod3.sed <- lmer(lsediment_kg_m3 ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod3.sed)
vif(sed.hlm.mod3.sed)
rand(sed.hlm.mod3.sed)

anova(sed.hlm.mod2.sed, sed.hlm.mod3.sed)

sink('./results/sediment_PAs.txt')
print(summary(sed.hlm.mod3.sed))

sink('./results/sediment_PAs_VIF.txt')
print(vif(sed.hlm.mod3.sed))

sink()
write.csv(round(summary(sed.hlm.mod3.sed)$coefficients,3), './results/sediment_PAs_coef.csv')


# (3-2-1) sediment, non-PAs
sed.hlm.mod3.sed.NPA <- lmer(lsediment_kg_m3 ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod3.sed.NPA)
vif(sed.hlm.mod3.sed.NPA)
rand(sed.hlm.mod3.sed.NPA)

sink('./results/sediment_nonPAs.txt')
print(summary(sed.hlm.mod3.sed.NPA))
sink()

write.csv(round(summary(sed.hlm.mod3.sed.NPA)$coefficients,3), './results/sediment_nonPAs_coef.csv')

# (3-2-2) sediment, PAs + Non_PAs
sed.hlm.mod3.sed.all <- lmer(lsediment_kg_m3 ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp  + (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod3.sed.all)
vif(sed.hlm.mod3.sed.all)
rand(sed.hlm.mod3.sed.all)

sink('./results/sediment_all.txt')
print(summary(sed.hlm.mod3.sed.all))
sink()

write.csv(round(summary(sed.hlm.mod3.sed.all)$coefficients,3), './results/sediment_all_coef.csv')

# (3-2-3) sediment, dam capacity instead of dam density
sed.hlm.mod3.sed.DC <- lmer(lsediment_kg_m3 ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program +lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), sed.hlm.sel)
summary(sed.hlm.mod3.sed.DC)
vif(sed.hlm.mod3.sed.DC)
rand(sed.hlm.mod3.sed.DC)

sink('./results/sediment_PAs_dam_capacity.txt')
print(summary(sed.hlm.mod3.sed.DC))
sink()

write.csv(round(summary(sed.hlm.mod3.sed.DC)$coefficients,3), './results/sediment_PAs_dam_capacity_coef.csv')
