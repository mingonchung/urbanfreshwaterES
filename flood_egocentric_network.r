rm(list=ls())

library(lme4)
library(lmerTest)
library(car)
library(memisc)
library(MuMIn)

setwd("D:/Min_Gon_Chung/Graduate_PhD/Research/Freshwater_data_030119/7_Data_analysis/process_040921/hlm")

flood.hlm <- read.csv("flood_watersheds_vars_gmis.csv", header=T, stringsAsFactors=FALSE)

### urban population > 300k, 2000-2010
flood.hlm.sel <- flood.hlm[which(flood.hlm$urb_pop_00_10 > 300),]

###########################
### log transformation ####
###########################
### (DV1) flood
flood.hlm.sel$lflood_city_pct <- log1p(flood.hlm.sel$flood_city_pct)

### (IV1) % of Forest, PAs
flood.hlm.sel$lforest_pct_PA_basin <- log1p(flood.hlm.sel$forest_pct_PA_basin)

### (IV1) % of Forest, no PAs
flood.hlm.sel$lforest_pct_noPA_basin <- log1p(flood.hlm.sel$forest_pct_noPA_basin)
flood.hlm.sel$lforest_pct_noPA_basin[is.na(flood.hlm.sel$lforest_pct_noPA_basin)] <- 0

### (IV1) % of Forest, PAs + no PAs
flood.hlm.sel$lforest_cover_pct <- log1p(flood.hlm.sel$forest_cover_pct)

### (IV2) % of Wetland, PAs
flood.hlm.sel$lwetland_pct_PA_basin <- log1p(flood.hlm.sel$wetland_pct_PA_basin)

### (IV2) % of Wetland, no PAs
flood.hlm.sel$lwetland_pct_noPA_basin <- log1p(flood.hlm.sel$wetland_pct_noPA_basin)
flood.hlm.sel$lwetland_pct_noPA_basin[is.na(flood.hlm.sel$lwetland_pct_noPA_basin)] <- 0

### (IV2) % of Wetland, PAs + no PAs
flood.hlm.sel$lwetland_cover_pct <- log1p(flood.hlm.sel$wetland_cover_pct)

### (IV3) dam density
flood.hlm.sel$ldam_density <- log1p(flood.hlm.sel$dam_density)

### (IV3) dam capacity
flood.hlm.sel$ldam_capacity <- log1p(flood.hlm.sel$dam_capacity)

### (IV4) % of irrigation
flood.hlm.sel$lirrigation_pct <- log1p(flood.hlm.sel$irrigation_pct)

### (IV5) Watershed areas, km2
flood.hlm.sel$lbasin_area <- log1p(flood.hlm.sel$basin_area)

### (IV6) Watershed distance, km
flood.hlm.sel$ldistance <- log1p(flood.hlm.sel$distance)

### (IV7) Elevation, meter
#min(flood.hlm.sel$elevation)
flood.hlm.sel$elevation_0 <- flood.hlm.sel$elevation + 36.59751
flood.hlm.sel$lelevation <- log1p(flood.hlm.sel$elevation_0)

### (IV8) Slope, degree
flood.hlm.sel$lslope <- log1p(flood.hlm.sel$slope)

### (IV9) IWS, 0-1
flood.hlm.sel$lIWS_program <- log1p(flood.hlm.sel$IWS_program)

### (IV10) Urban population, 1000 persons
flood.hlm.sel$lurb_pop_00_10 <- log1p(flood.hlm.sel$urb_pop_00_10)

### (IV11) Urban GDP per grid, billion USD in 2005 per degree (50km by 50km)
flood.hlm.sel$lurb_GDP_00_10 <- log1p(flood.hlm.sel$urb_GDP_00_10)

### (IV12) Temperature, Celsius degree
#min(flood.hlm.sel$city_temp)
flood.hlm.sel$city_temp_0 <- flood.hlm.sel$city_temp + 1.771133
flood.hlm.sel$lcity_temp <- log1p(flood.hlm.sel$city_temp_0)

### (IV13) Precipitation, mm
flood.hlm.sel$lcity_prcp <- log1p(flood.hlm.sel$city_prcp)

### (IV14) Impervious surface, pct
flood.hlm.sel$lcity_imper_pct <- log1p(flood.hlm.sel$city_imper_pct)

##########################
## Model 1: Null model ###
##########################
# (1-1) flood
flood.hlm.mod1 <- lmer(lflood_city_pct ~ 1 + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod1)
rand(flood.hlm.mod1)


####################################
## Model 2: Watersheds predictors ##
####################################
# (2-4) flood, PAs
flood.hlm.mod2 <- lmer(lflood_city_pct ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod2)
vif(flood.hlm.mod2)
rand(flood.hlm.mod2)


# (2-4-1) flood, non-PAs
flood.hlm.mod2.NPA <- lmer(lflood_city_pct ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod2.NPA)
vif(flood.hlm.mod2.NPA)
rand(flood.hlm.mod2.NPA)

# (2-4-2) flood, PAs + Non_PAs
flood.hlm.mod2.all <- lmer(lflood_city_pct ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod2.all)
vif(flood.hlm.mod2.all)
rand(flood.hlm.mod2.all)

# (2-4-3) flood, dam capacity instead of dam density
flood.hlm.mod2.DC <- lmer(lflood_city_pct ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod2.DC)
vif(flood.hlm.mod2.DC)
rand(flood.hlm.mod2.DC)


#############################################
## Model 3: watersheds + urban predictors ###
#############################################
# (3-4) flood, PAs
flood.hlm.mod3 <- lmer(lflood_city_pct ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod3)
vif(flood.hlm.mod3)
rand(flood.hlm.mod3)

anova(flood.hlm.mod2, flood.hlm.mod3)

sink('./results/flood_PAs.txt')
print(summary(flood.hlm.mod3))

sink('./results/flood_PAs_VIF.txt')
print(vif(flood.hlm.mod3))

sink()
write.csv(round(summary(flood.hlm.mod3)$coefficients,3), './results/flood_PAs_coef.csv')


# (3-4-1) flood, non-PAs
flood.hlm.mod3.NPA <- lmer(lflood_city_pct ~ lforest_pct_noPA_basin + lwetland_pct_noPA_basin + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod3.NPA)
vif(flood.hlm.mod3.NPA)
rand(flood.hlm.mod3.NPA)

sink('./results/flood_nonPAs.txt')
print(summary(flood.hlm.mod3.NPA))
sink()

write.csv(round(summary(flood.hlm.mod3.NPA)$coefficients,3), './results/flood_nonPAs_coef.csv')

# (3-4-2) flood, PAs + Non_PAs
flood.hlm.mod3.all <- lmer(lflood_city_pct ~ lforest_cover_pct + lwetland_cover_pct + ldam_density + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp  + (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod3.all)
vif(flood.hlm.mod3.all)
rand(flood.hlm.mod3.all)

sink('./results/flood_all.txt')
print(summary(flood.hlm.mod3.all))
sink()

write.csv(round(summary(flood.hlm.mod3.all)$coefficients,3), './results/flood_all_coef.csv')

# (3-4-3) flood, dam capacity instead of dam density
flood.hlm.mod3.DC <- lmer(lflood_city_pct ~ lforest_pct_PA_basin + lwetland_pct_PA_basin + ldam_capacity + lirrigation_pct + lbasin_area  + ldistance + lelevation + lslope + lIWS_program + lcity_imper_pct + lurb_pop_00_10 + lurb_GDP_00_10 + lcity_temp + lcity_prcp +  (1|City_ID), flood.hlm.sel)
summary(flood.hlm.mod3.DC)
vif(flood.hlm.mod3.DC)
rand(flood.hlm.mod3.DC)

sink('./results/flood_PAs_dam_capacity.txt')
print(summary(flood.hlm.mod3.DC))
sink()

write.csv(round(summary(flood.hlm.mod3.DC)$coefficients,3), './results/flood_PAs_dam_capacity_coef.csv')
