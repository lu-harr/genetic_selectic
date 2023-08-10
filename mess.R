library(raster)
library(viridisLite)
library(scales)
library(dismo)
library(gdistance)

###############################################################################
# Setting myself up
#source("~/Desktop/knowlesi/pilot_sites/code/0_health_sites_master.R")
setwd("~/Desktop/knowlesi")
source("~/Desktop/knowlesi/pilot_sites/code/1_pilot_site_setup_update.R")
# indices of districts of interest
DISTRICT_IND = c(2,3,4,5,6,7,9)

health_sites = readRDS(file = "~/Desktop/knowlesi/pilot_sites/health_sites.rds")

# catchment summary functions:
source(paste0(code_path,"5_catchment_summary.R"))
# apply catchment summary to distance and time catches
dist_catch_summary = lapply(health_sites$dist_catch,
                            catchment_summary_stats, lulc_covs)
time_catch_summary = lapply(health_sites$time_catch,
                            catchment_summary_stats, lulc_covs)
library(data.table)
dist_catch_summary = as.data.frame(data.table::rbindlist(dist_catch_summary))
time_catch_summary = as.data.frame(data.table::rbindlist(time_catch_summary))

source(paste0(code_path,"6_1_site_names.R"))

###############################################################################
# So, I've got sites, summarised for two types of catchments, 
# with an "obj" surface










