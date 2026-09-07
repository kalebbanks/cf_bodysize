#####################
#Script: 01_body_size_analysis.R
#Purpose: test different variables to see what best explains crawfish frog body size
#Inputs: data/sites, data/rasters
#Outputs: 
#Author: Kaleb M. Banks
#Date: 9/5/2026
#####################

#IMPORTANT: Skip to line 130 if you don't want to download all the rasters and want the ready to model dataframe

#####packages:
#install.packages("pgirmess")
#install.packages("geosphere")
#install.packages("geodata")
#install.packages("terra")
#install.packages("sf")
#install.packages("AICcmodavg")
#install.packages("nlme")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("glmm.hp")
#install.packages("ape")
library(pgirmess)
library(geosphere)
library(geodata)
library(terra)
library(sf)
library(AICcmodavg)
library(nlme)
library(dplyr)
library(ggplot2)
library(glmm.hp)
library(ape)


#####Download rasters

#Days above 5 degrees
download.file(
  url = "https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res01/CRUTS32/Hist/lt2_CRUTS32_Hist_8110.tif",
  destfile = "data/rasters/days_above_5.tif",
  mode = "wb"
)

#Days above 10 degrees
download.file(
  url = "https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res01/CRUTS32/Hist/lt3_CRUTS32_Hist_8110.tif",
  destfile = "data/rasters/days_above_10.tif",
  mode = "wb"
)

#Download bioclim annual temperature (1), temp seasonality (4), Annual rainfall (12), Precip seasonality (15)
worldclim_vars <- worldclim_global(
  var = "bio",
  res = 2.5,
  path = "data/rasters",
)


#####Load rasters
days_above_5 <- rast("data/rasters/days_above_5.tif")
days_above_10 <- rast("data/rasters/days_above_10.tif")
annual_temp <- rast("data/rasters/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif")
temp_seasonality <- rast("data/rasters/climate/wc2.1_2.5m/wc2.1_2.5m_bio_4.tif")
annual_rainfall <- rast("data/rasters/climate/wc2.1_2.5m/wc2.1_2.5m_bio_12.tif")
rainfall_seasonality <- rast("data/rasters/climate/wc2.1_2.5m/wc2.1_2.5m_bio_13.tif")


#####Load and format morpho data
cf_morpho <- read.csv("data/cf_morpho.csv")
#clean 
cf_morpho$Sex <- trimws(cf_morpho$Sex)
cf_morpho$Sex <- factor(cf_morpho$Sex, levels = c("F","M"))


#####Extract values from rasters
bio_vars <- c(
  annual_temp,
  temp_seasonality,
  annual_rainfall,
  rainfall_seasonality
)

names(bio_vars) <- c(
  "annual_temp",
  "temp_seasonality",
  "annual_rainfall",
  "rainfall_seasonality"
)

cf_points <- vect(
  cf_morpho,
  geom = c("Longitude", "Latitude"),
  crs = "EPSG:4326"
)

bio_vars_extract <- extract(bio_vars, cf_points)
days_above_10_extract <- extract(days_above_10, cf_points)
days_above_5_extract <- extract(days_above_5, cf_points)

cf_morpho <- cbind(
  cf_morpho,
  bio_vars_extract[, -1],                  
  days_above_10 = days_above_10_extract[, 2],
  days_above_5 = days_above_5_extract[, 2]
)


#####Create site column for random effect in lme
#there is like a million columns in this data frame, sorry blame Owen
cf_morpho <- cf_morpho %>%
  group_by(Latitude, Longitude) %>%
  mutate(Site = paste0("Site_", cur_group_id())) %>%
  ungroup()

cf_morpho <- cf_morpho %>%
  group_by(County) %>%
  mutate(county_site = paste0("Site_", cur_group_id())) %>%
  ungroup()

#One thing i am wondering if it would be better to group sites by coordinates (as in specimens with the exact same coordinates are treated as a site) or to group them by counties (all CF specimens in a county are considered a site). Grouping them by sites would greatly reduce the number of sites by half. As of right now I did the whole modeling framework with using coordinates as sites
length(unique(cf_morpho$Site))
length(unique(cf_morpho$county_site)) 
nrow(cf_morpho)         

#write.csv(cf_morpho, "data/cf_morpho_clean.csv")

#####Load and format Dataframe for modeling, can skip if you ran the above code ^
cf_morpho <- read.csv("data/cf_morpho_clean.csv")
cf_morpho$Sex <- factor(trimws(cf_morpho$Sex), levels = c("F","M"))

#####Hypothesis for body size differences in Anurans across wide environmental gradients
#These are taken from (Valenzuela-Sánchez et al., 2015)

#Hypothesis 1: Heat balance: this hypothesis predicts that larger individuals would be favoured in cold environments due to their reduced surface/mass ratio and enhanced thermal inertia. Variable: Annual mean temperature (temp_ann). Predicted effect on body size with decreased latitude: Negative

#Hypothesis 2:Temperature-size rule: this hypothesis predicts that larger body sizes are associated with colder climates, since a negative relationship between ontogenetic temperature and size at maturity has been found in many ectotherms. Variable: Annual mean temperature (temp_ann). Predicted effect on body size with decreased latitude: Negative. 

#Hypothesis 3: Optimal body temperature: assuming a constant optimal body temperature in all individuals of a given species and that thermoregulation is a critical factor in cold but not in warm environments (hypothesis formulated for squamate reptiles), this hypothesis predicts smaller body sizes to be associated with colder environments, as the increased surface/mass ratio of smaller individuals permit more rapid heating and cooling, improving thermoregulatory capacity in cooler climates. Conversely, relaxing selective pressure on surface/mass ratio in warm environments permits individuals attains larger sizes, increasing other size-related benefits. Variable: Annual mean temperature (temp_ann). Predicted effect on body size with decreased latitude: Positive.

#Hypothesis 4: Starvation resistance: larger body sizes are adaptive to more seasonal environments where individuals spend a long time in inactivity (e.g. hibernating) because larger individuals have higher resistance to starvation (i.e. energy stores increase with size faster than metabolic rate). Variable:Temperature seasonality (temp_seas). Predicted effect on body size with decreased latitude: Negative.

#Hypothesis 5:Growing season length:  larger body sizes are expected in less seasonal environments where conditions and resources allow a longer period of growth. Variable: For this, i looked at Lannoo's camera behavior paper (Stiles et al., 2017), looks like there is hardly any feeding behavior above temperatures of 5 degrees C, and at 10 degrees C activity really ramps up. I made variables below that are the days above 10/5 degrees C in a year.Predicted effect on body size with decreased latitude: Negative.

#Hypothesis 6: Water Availability: body sizes in amphibians are adaptive to drier environments because a lower surface/mass ratio reduces the loss of water. Variable:Temperature seasonality (m_rain_ann). Predicted effect on body size with decreased latitude: Negative.

#Hypothesis 7:Converse water availability:  larger body sizes in amphibians are associated with wetter climates because activity in this group is associated with high water availability (m_rain_ann). Predicted effect on body size with decreased latitude: Positive.

#These last two hypothesis aren't in Valenzunia-Sanchez but ones I wanted to test
#Hypothesis 8: Drought resistance Hypothesis: Large body sizes are more likely in area more prone to drought to better retain water. Variable: Precipation seasonality (m_rain_seas). Predicted effect on body size with decreased latitude: Negative.

#Hypothesis 9: Crawfish frogs are bigger in more northern areas because crayfish burrows are larger. Variable: Unfornuately there is no good way to test this. 

#Included latitude just for fun

#Papers mentioned: 
##Stiles, R. M., Halliday, T. R., Engbrecht, N. J., Swan, J. W., & Lannoo, M. J. (2017). Wildlife cameras reveal high resolution activity patterns in threatened crawfish frogs (Lithobates areolatus). Herpetological Conservation and Biology, 12, 160–170.
#Valenzuela-Sánchez, A., Cunningham, A. A., & Soto-Azat, C. (2015). Geographic body size variation in ectotherms: Effects of seasonality on an anuran from the southern temperate forest. Frontiers in Zoology, 12(1), 37. https://doi.org/10.1186/s12983-015-0132-y

#Look at correlations, expect most to be correlated with latitude
cor(cf_morpho[, c("Latitude","annual_temp","temp_seasonality",
                  "annual_rainfall","rainfall_seasonality",
                  "days_above_10","days_above_5")], use = "complete.obs")


# Linear Mixed-Effects models
m_null      <- lme(SVL_Ave ~ Sex,                           random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_lat       <- lme(SVL_Ave ~ Latitude + Sex,              random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_temp_ann  <- lme(SVL_Ave ~ annual_temp + Sex,           random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_temp_seas <- lme(SVL_Ave ~ temp_seasonality + Sex,      random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_rain_ann  <- lme(SVL_Ave ~ annual_rainfall + Sex,        random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_rain_seas <- lme(SVL_Ave ~ rainfall_seasonality + Sex,  random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_days10    <- lme(SVL_Ave ~ days_above_10 + Sex,         random = ~ 1 | Site, data = cf_morpho, method = "ML")
m_days5    <- lme(SVL_Ave ~ days_above_5 + Sex,         random = ~ 1 | Site, data = cf_morpho, method = "ML")

# List for AICc ranking
models <- list(
  "Null"                 = m_null,
  "Latitude"             = m_lat,
  "Annual Temp"          = m_temp_ann,
  "Temp Seasonality"     = m_temp_seas,
  "Annual Rainfall"      = m_rain_ann,
  "Rainfall Seasonality" = m_rain_seas,
  "Days Above 10"        = m_days10, 
  "Days Above 5"        = m_days5
  
)

# AICc Table
aictab(cand.set = models)

top_model_reml <- lme(
  SVL_Ave ~ temp_seasonality + Sex,
  random = ~ 1 | Site,
  data = cf_morpho,
  method = "REML"
)

summary(top_model_reml)
intervals(top_model_reml)




#####Moran's I test

# Moran's I tests whether values that are geographically close together are more similar than values that are far apart (spatial autocorrelation), more than expected by random chance. It returns one number (roughly -1 to +1): positive = nearby sites are similar (spatial clustering), near 0 = no spatial pattern, negative = nearby sites are more different than expected.

# Why we're running this: our top model explains SVL using temp_seasonality, a climate variable that is itself spatially structured (nearby sites share similar climate). If our model's leftover errors (residuals) still show a spatial pattern after accounting for temp_seasonality, that would suggest some other spatially-patterned variable we didn't measure (e.g. elevation, genetic structure, an unmodeled climate axis) is still influencing body size. A non-significant result on the residuals tells us the model has adequately captured the spatial signal in the data, and our p-values aren't being artificially inflated by unmodeled spatial autocorrelation.

# Moran's I test: Raw SVL
#We expect this to be significant, it just confirms there is spatial climate gradient in the data, we dont really care about this test. 
site_svl <- cf_morpho %>% distinct(Site, Latitude, Longitude, SVL_Ave) %>% arrange(Site)
dist_matrix <- distm(site_svl[, c("Longitude","Latitude")], fun = distHaversine)
inv_dist <- 1 / dist_matrix; diag(inv_dist) <- 0
diag(inv_dist) <- 0
#For sites that share the same coordinates I made their distances 1. This analysis needs values larger than 0. 
inv_dist[is.infinite(inv_dist)] <- 1

Moran.I(site_svl$SVL_Ave, inv_dist) 


#Moran's I test: Residuals of top model
#We expect this to come back non-significant, which would mean there is no leftover spatially patterned variable we're missing
cf_morpho$resid <- residuals(top_model_reml, type = "normalized")

Moran.I(cf_morpho$resid, inv_dist) 


#Raw SVL showed significant spatial autocorrelation (I = 0.363, p < 0.001); model residuals did not (I = 0.074, p = 0.228), indicating no residual spatial structure left unaccounted for by the top model.


#####Hierarchical partitioning

#It calculates the indepent explanatory power of each variable by averaging goodness of fit measures across all possible subset models. It solves problems caused by multicollinearity showing which variables truly explain the most variance independently. 

m_full_hp <- lme(SVL_Ave ~ annual_temp + temp_seasonality + days_above_5 +
                   annual_rainfall + rainfall_seasonality + Sex,
                 random = ~ 1 | Site,
                 data = cf_morpho, method = "ML")

hp_result <- glmm.hp(m_full_hp)



#####Figures

#Hierarchical partitioning figure
hp_df <- data.frame(
  variable = c("Annual Temp", "Temp Seasonality", "Days Above 5", "Sex", "rain seasonality", "annual rainfall"),
  I_perc = c(22.70, 31.06, 23.78, 5.60, 5.16, 11.70)
)
hp_df$variable <- factor(hp_df$variable, levels = hp_df$variable[order(-hp_df$I_perc)])
ggplot(hp_df, aes(x = variable, y = I_perc)) +
  geom_col(fill = "grey40", width = 0.6) +
  labs(x = NULL, y = "Independent effect (%)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


#Morin's I for mean adult body size and residuals
site_sf <- st_as_sf(site_svl, coords = c("Longitude", "Latitude"), crs = 4326)
site_albers <- st_transform(site_sf, crs = 5070)  # CONUS Albers Equal Area
coords_km <- st_coordinates(site_albers) / 1000  # meters to km

correlog_svl <- pgirmess::correlog(coords_km, site_svl$SVL_Ave, method = "Moran", nbclass = 10)
df_svl <- as.data.frame(correlog_svl)
colnames(df_svl) <- c("dist.class", "Moran.I", "p.value", "n")
df_svl$type <- "Body size"

cf_morpho_sf <- st_as_sf(cf_morpho, coords = c("Longitude", "Latitude"), crs = 4326)
cf_morpho_albers <- st_transform(cf_morpho_sf, crs = 5070)
coords_km_ind <- st_coordinates(cf_morpho_albers) / 1000

correlog_resid <- pgirmess::correlog(coords_km_ind, cf_morpho$resid, method = "Moran", nbclass = 10)
df_res <- as.data.frame(correlog_resid)
colnames(df_res) <- c("dist.class", "Moran.I", "p.value", "n")
df_res$type <- "Residuals"

correlog_combined <- rbind(df_res, df_svl)

ggplot(correlog_combined, aes(x = dist.class, y = Moran.I,
                              shape = type, fill = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(group = type), color = "black") +
  geom_point(size = 3, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c("Body size" = 21, "Residuals" = 21)) +
  scale_fill_manual(values = c("Body size" = "black", "Residuals" = "white")) +
  labs(x = "Distance class (km)", y = "Moran's I", shape = NULL, fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")


#SVL with Temperature seasonality using our most supported model
newdata <- data.frame(
  temp_seasonality = seq(min(cf_morpho$temp_seasonality, na.rm = TRUE),
                         max(cf_morpho$temp_seasonality, na.rm = TRUE), length.out = 100),
  Sex = factor(levels(cf_morpho$Sex)[1], levels = levels(cf_morpho$Sex))
)
newdata$pred <- predict(top_model_reml, newdata = newdata, level = 0)

ggplot() +
  geom_point(data = cf_morpho,
             aes(x = temp_seasonality, y = SVL_Ave, shape = Sex, fill = Subspecies),
             size = 5, color = "black", alpha = 1) +
  geom_line(data = newdata, aes(x = temp_seasonality, y = pred),
            color = "firebrick", linewidth = 1) +
  scale_shape_manual(values = c("M" = 21, "F" = 24)) +
  scale_fill_manual(values = c("Intergrade" = "orange","Southern" = "forestgreen", "Northern" = "steelblue")) +
  labs(x = "Temperature seasonality", y = "Snout-vent length (mm)",
       shape = "Sex", fill = "Subspecies") +
  theme_classic(base_size = 13)

