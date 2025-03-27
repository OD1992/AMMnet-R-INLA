rm(list=ls())
###run this for package installation
packages <- c("malariaAtlas", "raster", "sp", "tidyverse",
              "lattice", "gridExtra", "devtools", "rlang")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) { 
  install.packages(setdiff(packages, rownames(installed.packages()))) }

#For INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

##load packagaes
library(INLA)
library(malariaAtlas)
library(raster)
library(sp)
library(tidyverse)
library(lattice)     
library(gridExtra)

# Data Cleaning -----------------------------------------------------------

#get data
MDG_shp <- getShp(ISO = "MDG", admin_level = c("admin0"))
prj <- crs(MDG_shp)
e <- extent(MDG_shp)
# MDG_pr_data <- getPR(country = 'Madagascar', species = 'both')
# autoplot(MDG_pr_data, facet = FALSE)
#pull the new Madagascar data
MDG_pr_data <- read_csv("/Volumes/map_per/GBD2022/Processing/Stages/08_PR_DB_Import_Export/Checkpoint_Outputs/20220908/pfpr.csv") %>%
  filter(country == "Madagascar") %>% 
  filter(year_start > 2020) %>% 
  filter(!is.na(pf_pos) & !is.na(examined))

y <- round(MDG_pr_data$pf_pos)
n <- MDG_pr_data$examined
# Covariates cleaning -----------------------------------------------------
#This section will load the covarites; clip them to your analysis shapefile and scale them
#this is then outputted to your covariates folder
#NOTE: the raw covariates are available from MAP upon request; else you will have to modify this
# for your own covariate rasters


 # temp15 <- raster('/Volumes/map_per/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Annual/LST_Day_v6.2021.Annual.mean.5km.mean.tif') %>%  # temp15 <- raster('/Volumes/map_per/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Annual/LST_Day_v6.2021.Annual.mean.5km.mean.tif') %>% 
 #   crop(MDG_shp) %>% mask(MDG_shp) %>% scale() %>% 
 #   writeRaster( file = "covariates/LST_day.tif", overwrite=T)
 # plot(temp15)
 # 
 # access_hf <- raster('/Volumes/map_per/mastergrids/Other_Global_Covariates/Accessibility/Weiss/2019/accessibility_to_healthcare_2019.tif') %>% 
 #   aggregate(fact=5)%>% crop(MDG_shp) %>% mask(MDG_shp) %>%  scale()%>% 
 #   writeRaster( file = "covariates/Access_hf.tif", overwrite=T)
 # plot(access_hf)
 # 
 # rain <- raster('/Volumes/map_per/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/5km/Annual/chirps-v2-0.2021.Annual.sum.5km.NN.tif')%>% 
 #   crop(MDG_shp) %>% mask(MDG_shp) %>% scale()%>% 
 #   writeRaster( file = "covariates/Rain.tif", overwrite=T)
 # plot(rain)
 # 
 # evi <- raster('/Volumes/map_per/mastergrids/Other_Global_Covariates/Aridity_v2/5km/Synoptic/Aridity_Index_v2.Synoptic.Overall.Data.5km.mean.tif')%>% 
 #   crop(MDG_shp) %>% mask(MDG_shp) %>% scale()%>% 
 #   writeRaster( file = "covariates/EVI.tif", overwrite=T)
 # plot(evi)
 # 
 # elevation <- raster('/Volumes/map_per/mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif')%>% 
 #   crop(MDG_shp) %>% mask(MDG_shp) %>% scale()%>% 
 #   writeRaster( file = "covariates/Elevation.tif", overwrite=T)   # writeRaster( file = "covariates/Elevation.tif", overwrite=T)
 # plot(elevation)


# Covariates cleaning II -----------------------------------------------------
covs.list <- list.files('covariates', pattern = ".tif$", full.names = T)
covs <- stack(covs.list)
plot(covs)

#extract the covariates
#first create spatial points dataframe to extract
MDG_points <- MDG_pr_data[,c("longitude", "latitude")] %>% 
  as.matrix %>% 
  SpatialPoints(proj4string=prj)

#use the extract function from raster to extract the covariates for each point
covs_df <- raster::extract(covs, MDG_points)

#combine covariates with data
MDG_pr_data = cbind(MDG_pr_data, covs_df)
#View(MDG_pr_data) #look at the data

# Data Export ----------------------------------------------------
#write out cleaned data
save(y, n, cov, covs_df, covs, MDG_pr_data, MDG_shp, file = 'input/MDG_clean.Rdata')
