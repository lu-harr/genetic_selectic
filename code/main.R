# working from jev_multi_site
getwd()
setwd("~/Desktop/knowlesi/multi_site")

library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(terra)
library(RColorBrewer)
library(scales)
library(viridisLite)
library(data.table)
library(rPref)
library(idpalette)

AGG_FACTOR <- 10

#source("code/iterative_select_funcs.R")
pinks <- colorRampPalette(colors = c("#f7f7f7", "#c23375"))
purps <- colorRampPalette(brewer.pal(9,"Purples"))
#extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(100)

greens = colorRampPalette(colors = c("#e0fbb7","#8ACE00","#6DCD59FF","#2f9a16"))
apple = "#6DCD59FF"
brat = "#8ACE00"
berry = "#c23375"
purp = "#8612ff"
oranges = colorRampPalette(brewer.pal(9,"Oranges"))

subfigure_label = function(plot_region, x_displacement, y_displacement, label,
                           cex.label=1){
  # gives me an (a)
  text(plot_region[1] + (plot_region[2] - plot_region[1])*x_displacement, 
       plot_region[3] + (plot_region[4] - plot_region[3])*y_displacement,
       label, cex=cex.label)
}

empty_plot_for_legend = function(){
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
}

states <- st_read("~/Desktop/jev/data/admin/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")
# all_mozzies <- read.csv('~/Desktop/jev/from_Freya_local/JEV_secure/data/national_mozzie_data/mosquito_detections_all_w_qld_sa.csv')

# check this is the data I end up using in JEV 1?
# all_mozzies <- read.csv('~/Desktop/jev/from_Freya_local/JEV_secure/data/phase_2/clean/pathogen/all_moz_surveillance.csv') %>%
#   rename(host_type=host_species) %>%
#   dplyr::select(
#     -c(site,
#        n_individual_spp_tests,
#     )
#   ) %>%
#   drop_na(longitude:virus_detection)
# all_mozzies$data_source[all_mozzies$data_source == "VIC"] = "Vic"

# potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
# 
# #potential_continuous <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
# hpop <- raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif')
# 
# vic_shp <- states %>%
#   filter(STE_NAME21 == "Victoria") %>%
#   st_simplify(dTolerance = 1000)
# 
# vic_shadow <- vic_shp %>%
#   st_buffer(dist = 5000) %>% # 80% certain we're in metres
#   st_simplify(dTolerance = 5000)
# # ugly ! ah well !
# 
# potential_vic_buffered <- potential_continuous %>%
#   crop(vic_shadow) %>%
#   raster::mask(vic_shadow)
# # probably a good idea to check northern catches only pick up pixels in buffer
# 
# potential_vic <- potential_continuous %>%
#   crop(vic_shadow) %>%
#   raster::mask(vic_shp)
# 
# hpop_vic_buff <- hpop %>%
#   crop(vic_shadow) %>%
#   raster::mask(vic_shp)
# 
# AGG_FACTOR <- 10
# 
# potential_guelphia <- potential_continuous %>%
#   crop(guelphia_extent) %>%
#   aggregate(AGG_FACTOR) %>%
#   t()
# 
# guelphia_hpop <- hpop %>%
#   crop(guelphia_extent) %>%
#   aggregate(AGG_FACTOR) %>%
#   t()
# 
# wtmat = matrix(1, nrow=3, ncol=3)
# blurred_potential = focal(potential_guelphia, wtmat, sum, na.rm=TRUE)
# 
# site_ids = as.data.frame(rasterToPoints(blurred_potential))
# site_ids$id = which(!is.na(values(blurred_potential)))
# nrow(site_ids)
# 
# single_site_quants = quantile(potential_guelphia, probs=seq(0,1,0.05))
# sel = potential_guelphia
# sel[sel < single_site_quants[19]] = NA
# 
# plot(potential_guelphia)
# plot(sel, col="red", add=TRUE)
# 
# distance_matrix = rasterToPoints(blurred_potential)[,c("x","y")] %>%
#   dist(upper=TRUE) %>%
#   as.matrix()
# 
# obj_stack <- stack(potential_guelphia, 
#                    guelphia_hpop,
#                    blurred_potential)
# names(obj_stack) <- c("potential", "hpop", "catch_potential")
# site_ids <- obj_stack %>%
#   mask(blurred_potential) %>%
#   rasterToPoints() %>%
#   as.data.frame()
# site_ids$id <- which(!is.na(values(blurred_potential)))
# 
# # objective_list = list(
# #   obj1 = function(site_ids, 
# #                   dist_mat = matrix(NA), 
# #                   raster_stack = stack(), 
# #                   pix_weights = c()){
# #     # sum objective values of pix in catch (first stack in raster)
# #     sum(raster_stack[["catch_potential"]][site_ids])
# #   },
# #   obj2 = function(site_ids, 
# #                   dist_mat = matrix(NA), 
# #                   raster_stack = stack(), 
# #                   pix_weights = c()){
# #     # sum objective values of pix in catch (second stack in raster)
# #     sum(raster_stack[["hpop"]][site_ids])
# #   },
# #   # this one is supposed to be network distance (Euclidean, MST via Prim)
# #   obj3 = function(site_ids = c(),
# #                   dist_mat = matrix(NA),
# #                   raster_stack = stack(),
# #                   pix_weights = c()){
# #     # this one only needs site_ids and dist_mat!
# #     if (length(site_ids) == 2){
# #       return(dist_mat[site_ids[1], site_ids[2]])
# #     }
# #     site_ids = sort(unlist(site_ids)) # does this need to happen?
# #     dist_mat = dist_mat[site_ids, site_ids]
# #     tmp_mst = ape::mst(as.dist(dist_mat))
# #     picks = cbind(expand.grid(rownames(tmp_mst), colnames(tmp_mst)), 
# #                   as.vector(tmp_mst))
# #     picks[, 1:2] = cbind(as.numeric(picks[, 1]), as.numeric(picks[, 2]))
# #     picks = picks[picks[, 3] == 1,]
# #     picks[, 1:2] = t(apply(picks[, 1:2], 1, sort))
# #     picks = unique(picks)
# #     sum(dist_mat[as.matrix(picks[, 1:2])])
# #   }
# # )
# 
# 
# 
# 
# 
# 
# 
# 
# 
