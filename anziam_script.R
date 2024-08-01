# theoretically I *could* have a go at JEV mozzies
# want to maximise network distance .. 
# could define network in grid which would get me manhattan distance ..
# could give sites a radius and minimise overlap - this might be the MVP

# I think going back to the 10kmsq version might be less hectic ... optimise distance ...?
# need to fix mozzie/bird masks to remove IDN/PNG

library(sf)
library(raster)
library(dplyr)
library(terra)
library(RColorBrewer)
library(scales)
source("~/Desktop/knowlesi/multi_site/iterative_select_funcs.R")

path <- "~/Desktop/jev/from_Freya_local/JEV_secure/"
# need mozzies for current trapping
all_mozzies <- read.csv('~/Desktop/jev/from_Freya_local/JEV_secure/data/national_mozzie_data/mosquito_detections_all_w_qld_sa.csv')

AGG_FACTOR <- 10

potential_continuous <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
#potential_continuous <- aggregate(potential_continuous, fact=10, fun='mean')
#overall_power <- raster('output/national_detection_probability_cut.tif')
states = st_read("~/Desktop/jev/data/admin/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")

guelphia_extent <- c(141, 146, -37.9,-32.9)

# vic_shp <- states %>%
#   filter(STE_NAME21 == "Victoria") %>%
#   st_simplify(dTolerance = 1000)

# vic_shadow <- vic_shp %>%
#   st_buffer(dist = 5000) %>% # 80% certain we're in metres
#   st_simplify(dTolerance = 5000)
# ugly ! ah well !

# potential_vic_buffered <- potential_continuous %>%
#   crop(vic_shadow) %>%
#   raster::mask(vic_shadow)
# probably a good idea to check northern catches only pick up pixels in buffer

# potential_vic <- potential_continuous %>%
#   crop(vic_shp) %>% # was shadow
#   raster::mask(vic_shp)

potential_guelphia <- potential_continuous %>%
  crop(guelphia_extent) %>%
  aggregate(AGG_FACTOR) %>%
  t()

plot(potential_continuous)
# oops this is rotated now :)
plot(t(potential_guelphia), col="red", add=TRUE)
plot(potential_guelphia)

# vic_hpop <- raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif') %>%
#   crop(vic_shp) %>% # was shadow
#   raster::mask(vic_shp) %>%
#   aggregate(AGG_FACTOR)

# Perhaps this gives me away? Could move the "city" northwards?
guelphia_hpop <- raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif') %>%
  crop(guelphia_extent) %>%
  aggregate(AGG_FACTOR) %>%
  t()

# Done in buffered space to pick up land pixels in SA/NSW (mosquitoes don't see borders)
# wtmat = focalWeight(guelphia_potential, 0.1, type="circle")
# wtmat[wtmat != 0] = 1
wtmat = matrix(1, nrow=3, ncol=3)
blurred_potential = focal(potential_guelphia, wtmat, sum, na.rm=TRUE)
# helpfully omits border

# need id in buffered space, but only want victorian pixels
site_ids = as.data.frame(rasterToPoints(blurred_potential))
site_ids$id = which(!is.na(values(blurred_potential)))
nrow(site_ids)
# after aggregating, we're down from 327,652
# also, won't need to worry about catchments
# can throw away "buffer zone"

# good idea to check this :)
which(is.na(site_ids$continuous.suit.vectors.and.avian))

# single_site_quants = quantile(blurred_potential_masked, probs=seq(0,1,0.05))
single_site_quants = quantile(potential_guelphia, probs=seq(0,1,0.05))
sel = potential_guelphia
sel[sel < single_site_quants[19]] = NA

{png("~/Desktop/jev/anziam24/assess_guelphia.png",
    height=900,
    width=2000,
    pointsize=30)
par(bty="n", xpd=NA, mar=c(2.1,2.1,2.1,2.1), mfrow=c(1,2))
pn_cols <- colorRampPalette(colors = c("#f7f7f7", "#c23375"))(1000)
purps <- colorRampPalette(brewer.pal(9,"Purples"))(1000)
# plot(potential_vic, col=brewer.pal(9, "Greys"), legend=FALSE,
#      axes=FALSE)#, xlim=c(140.9,141.2), ylim=c(-34.15,-33.9))
plot(potential_guelphia, col=pn_cols, axes=FALSE)
legend(x=145, y=-33.5, pch=1, "Existing surveillance")
legend(x=145, y=-34.5, fill="orange", "Area at risk")
par(new=TRUE)
# plot(potential_vic, col=pn_cols, axes=FALSE)#, add=TRUE)#, xlim=c(140.9,141.2), ylim=c(-34.15,-33.9))
# plot(st_geometry(vic_shp), add=TRUE)
par(new=TRUE)
plot(sel, col=alpha("orange", 0.6), new=TRUE, legend=FALSE, axes=FALSE)

plot(sqrt(guelphia_hpop), col=purps, legend=FALSE, axes=FALSE)
# plot(st_geometry(vic_shp), add=TRUE)
par(new=TRUE)
# this does look suspiciously like the Murray .. maybe I just wnat to pick somewhere in the Kimberly ?
plot(sel, col=alpha("orange", 0.6), new=TRUE, legend=FALSE, axes=FALSE)
# points(all_mozzies[all_mozzies$state == "Vic",c("longitude", "latitude")], cex=0.6)
dev.off()}


#########################################################################
# OBJECTIVE LIST WITH DISTANCE MATRIX

distance_matrix = rasterToPoints(blurred_potential)[,c("x","y")] %>%
  dist(upper=TRUE) %>%
  as.matrix()

dim(distance_matrix)

plot(raster(distance_matrix))

obj_stack <- stack(potential_guelphia, 
                   guelphia_hpop,
                   blurred_potential)
names(obj_stack) <- c("potential", "hpop", "catch_potential")
site_ids <- obj_stack %>%
  mask(blurred_potential) %>%
  rasterToPoints() %>%
  as.data.frame()
site_ids$id <- which(!is.na(values(blurred_potential)))

dim(site_ids)

# Here's an example objective_list
# assuming no NA pixels in study area !
objective_list = list(
  obj1 = function(site_ids, 
                  dist_mat = matrix(NA), 
                  raster_stack = stack(), 
                  pix_weights = c()){
    # sum objective values of pix in catch (first stack in raster)
    sum(raster_stack[["catch_potential"]][site_ids])
  },
  obj2 = function(site_ids, 
                  dist_mat = matrix(NA), 
                  raster_stack = stack(), 
                  pix_weights = c()){
    # sum objective values of pix in catch (first stack in raster)
    sum(raster_stack[["hpop"]][site_ids])
  },
  # this one is supposed to be network distance (Euclidean, MST via Prim)
  obj3 = function(site_ids = c(),
                  dist_mat = matrix(NA),
                  raster_stack = stack(),
                  pix_weights = c()){
    # this one only needs site_ids and dist_mat!
    if (length(site_ids) == 2){
      return(dist_mat[site_ids[1], site_ids[2]])
    }
    site_ids = sort(unlist(site_ids)) # does this need to happen?
    dist_mat = dist_mat[site_ids, site_ids]
    tmp_mst = ape::mst(as.dist(dist_mat))
    picks = cbind(expand.grid(rownames(tmp_mst), colnames(tmp_mst)), 
                  as.vector(tmp_mst))
    picks[, 1:2] = cbind(as.numeric(picks[, 1]), as.numeric(picks[, 2]))
    picks = picks[picks[, 3] == 1,]
    picks[, 1:2] = t(apply(picks[, 1:2], 1, sort))
    picks = unique(picks)
    sum(dist_mat[as.matrix(picks[, 1:2])])
  }
)

non_raster_objective_list = list(
  obj1 = function(site_ids, 
                  dist_mat = matrix(NA), 
                  vector_to_be_indexed = c()){
    # sum objective values of pix in catch (first stack in raster)
    sum(vector_to_be_indexed[site_ids], na.rm=TRUE)
  },
  # this one is supposed to be network distance (Euclidean, MST via Prim)
  obj2 = function(site_ids = c(),
                  dist_mat = matrix(NA),
                  vector_to_be_indexed = c()){
    # this one only needs site_ids and dist_mat!
    if (length(site_ids) == 2){
      return(dist_mat[site_ids[1], site_ids[2]])
    }
    site_ids = sort(unlist(site_ids)) # does this need to happen?
    dist_mat = dist_mat[site_ids, site_ids]
    tmp_mst = ape::mst(as.dist(dist_mat))
    picks = cbind(expand.grid(rownames(tmp_mst), colnames(tmp_mst)), 
                  as.vector(tmp_mst))
    picks[, 1:2] = cbind(as.numeric(picks[, 1]), as.numeric(picks[, 2]))
    picks = picks[picks[, 3] == 1,]
    picks[, 1:2] = t(apply(picks[, 1:2], 1, sort))
    picks = unique(picks)
    sum(dist_mat[as.matrix(picks[, 1:2])])
  }
)

#choose(nrow(site_ids), 2) # 6 mill ...

# this takes way too long for what it is ...

# {t1 = Sys.time()
# obj1_choose2 <- combn(100, 2, objective_list[[1]], TRUE, 
#                       raster_stack = obj_stack)
# t2 = Sys.time()}
# 
# {t1 = Sys.time()
#   obj2_choose2 <- combn(100, 2, objective_list[[2]], TRUE, 
#                         raster_stack = obj_stack)
#   t2 = Sys.time()}
# 
# (t2-t1)/choose(100, 2)*choose(3600, 2)/60
# 
# choose(1000, 2)
# 
# {t1 = Sys.time()
# obj3_choose2 <- combn(100, 2, objective_list[[3]], TRUE,
#                       dist_mat = distance_matrix)
# t2 = Sys.time()}
# 
# (t2-t1)/choose(100, 2)*choose(3600, 10)/60

# but wait! Objectives 1 and 2 do not require simultaneous selection of sites: sites can be selected successively !
# (assuming mozzies can't fly between pixels)

# So let's do some Pareto-optimisation of objectives 1 and 2 !
# (Can probably throw out things that aren't super great under objectives 1 and 2)

# Approach 1: look exhaustively at best n sites via sapply and sample?
best_sites <- order(site_ids$potential, decreasing = TRUE)[1:100] # removed indexing here
head(best_sites)
plot(obj_stack$potential, main="'Guelphia' potential transmission suitability")
text(site_ids[best_sites[1:22], c("x","y")], labels=1:22, col="blue", cex=0.8)
plot(sqrt(obj_stack$hpop), main="'Guelphia' hpop")

# {t1 = Sys.time()
# tmp2 <- sapply(1:100, function(x){
#   sum(site_ids$potential[sample(best_sites$id, 10, replace=TRUE)])
# })
# t2 = Sys.time()}
# (t2-t1)/100*choose(50, 8)/60

id_ras <- trim(blurred_potential)
#id_ras[!is.na(id_ras)] <- which(!is.na(values(blurred_potential))) # case for NAs
id_ras[] <- 1:ncell(id_ras)
neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
# don't actually want *weights* per se .. but that might be useful later ...
neigh_mat[!neigh_mat == 0] = 1
# jumping into terra because of option for fun argument that returns vector
# this gives me a layered SpatRaster
catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c) # so much quicker than my diy version
# need to remove zero layers ...
catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
# and remove non-victorian pixels that were part of the buffer ...
catchment_stack <- mask(catchment_stack, rast(trim(blurred_potential)))
# now I want this as a list of vectors? Or a matrix?
catch_membership_mat <- values(catchment_stack, mat=TRUE)
# this leaves NAs in there ... better remove those
catch_membership_mat <- catch_membership_mat[!is.na(values(trim(blurred_potential))),]

# Approach 2: sample designs from best n sites
# TIL combn() rejects 100 choose 10
{t1 = Sys.time()
  tmp = combn(22,10)
  #tmp = tmp[, sample.int(ncol(tmp), 100)]
  t2 = Sys.time()
  obj1_choose10 <- apply(tmp, 2, function(x){
    #c(x, sum(site_ids$catch_potential[best_sites[x]])) # removed "id" here
    # can't count sites in more than one catchment twice:
    c(x, sum(site_ids$potential[unique(as.vector(catch_membership_mat[best_sites[x],]))], na.rm=TRUE))
  })
  t3 = Sys.time()
}

(t3-t1)/choose(22,10)*choose(3600,10)/60/60/24/365 # time to do this for all sites in years
obj1_choose10 = t(obj1_choose10)
dim(obj1_choose10)
plot(obj1_choose10[, 11]) # objective
saveRDS(obj1_choose10, "~/Desktop/jev/anziam24/obj1_choose10.rds")
obj1_choose10 <- readRDS("~/Desktop/jev/anziam24/obj1_choose10.rds")

# Approach 3: go straight for combn?
# {t1 = Sys.time()
# obj1_choose10 <- combn(100, 2, function(x){sum(site_ids$potential[best_sites[x, "id"]])}, 
#                        TRUE)
# t2 = Sys.time()}
# (t2-t1)/choose(100, 2)*choose(100, 10)/60

# objective 2
# which sites to choose? well I know what the minimum will be and I can work out roughly what the max will be ...

# lower bound (best case):
best_design_obj2_eg <- c(nrow(trim(blurred_potential)) + 1:3, 
                          nrow(trim(blurred_potential))*2 + 1:4, 
                          nrow(trim(blurred_potential))*3 + 1:3)
plot(blurred_potential, ylim=c(145.62, 145.8), xlim=c(-37.78, -37.52))
points(site_ids[best_design_obj2_eg, c("x","y")])
objective_list[[3]](site_ids = best_design_obj2_eg, dist_mat = distance_matrix)

# upper bound:
# look at https://dickbrus.github.io/SpatialSamplingwithR/RegularGridSpatialCoverage.html
# this is getting silly ... I'll just overlay a grid
worst_design_ish_obj2_eg <- sampleRegular(trim(blurred_potential), size=10, xy=TRUE, cells=TRUE)
objective_list[[3]](site_ids = worst_design_ish_obj2_eg[,"cell"], dist_mat = distance_matrix)

# compare best designs under first objective to second objective?
# {t1 = Sys.time()
#   tmp = combn(13,10) # check this!
#   #tmp = tmp[, sample.int(ncol(tmp), 100)]
#   t2 = Sys.time()
#   obj2_bestobj1_choose10 <- apply(tmp, 2, function(x){
#     c(x, objective_list[[3]](site_ids = best_sites[x], dist_mat = distance_matrix)) # removed "id" here
#   })
#   t3 = Sys.time()}
# saveRDS(obj2_bestobj1_choose10, "~/Desktop/jev/anziam24/obj2_bestobj1_choose10.rds")
obj2_bestobj1_choose10 <- readRDS("~/Desktop/jev/anziam24/obj2_bestobj1_choose10.rds")
# (t3-t1)/choose(13,10)*choose(3600, 10)/60/60/24/365
 
obj1_choose10 = cbind(obj1_choose10, obj2_bestobj1_choose10[11,])
tmp = obj1_choose10
dim(tmp)
# need to order by obj1 ..
tmp = tmp[order(tmp[,11], decreasing=TRUE),]
matplot(tmp[1:100, 1:10], type="l") # I'm a bit perturbed by this

pareto10 <- data.frame(obj1_choose10)
names(pareto10) <- c(paste("site", 1:10), "sum_risk", "net_dist")
nrow(pareto10)
#plot(pareto10$sum_risk, pareto10$net_dist)
pareto10 <- psel(pareto10, high("sum_risk")*low("net_dist"))
pareto10 <- pareto10[order(pareto10$sum_risk),]
#lines(pareto10$sum_risk, pareto10$net_dist, col="red")#, pch=16)


#obj1_choose10 = cbind(obj1_choose10, obj2_bestobj1_choose10[11,])
plot(obj1_choose10[,11], obj1_choose10[,12], 
     xlim=c(0, max(obj1_choose10[,11])),
     ylim=c(objective_list[[3]](site_ids = best_design_obj2_eg, dist_mat = distance_matrix),
            objective_list[[3]](site_ids = worst_design_ish_obj2_eg[,"cell"], dist_mat = distance_matrix)),
     xlab="Total potential risk",
     ylab="Network distance",
     cex.lab=1.3, cex=0.7,
     main="All designs of 20 individually riskiest sites", cex.main=1.4, type="n")
points(pareto10[,11], pareto10[,12], col="#6DCD59FF")
lines(pareto10[,11], pareto10[,12], col="#6DCD59FF", lwd=2)
# is pretty good no?
# yeah that looks right

# hmmm this didn't work how I thought it would... nope I found the bug false alarm


# I'm missing an opportunity for brat-themed figures here lmao


###############################################################################
# BELOW CODE DOES UN-AGGREGATED CASE:

id_ras <- potential_vic_buffered
id_ras[!is.na(id_ras)] <- which(!is.na(values(potential_vic_buffered)))
neigh_mat <- focalWeight(id_ras, 0.05, "circle")
# don't actually want *weights* per se .. but that might be useful later ...
neigh_mat[!neigh_mat == 0] = 1
# jumping into terra because of option for fun argument that returns vector
# this gives me a layered SpatRaster
catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c) # so much quicker than my diy version
# need to remove zero layers ...
catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
# and remove non-victorian pixels that were part of the buffer ...
catchment_stack <- mask(catchment_stack, rast(potential_vic))
# now I want this as a list of vectors? Or a matrix?
catch_membership_mat <- values(catchment_stack, mat=TRUE)
# this leaves NAs in there ... better remove those
catch_membership_mat <- catch_membership_mat[!is.na(values(potential_vic)),]

plot(raster(catch_membership_mat))

# objective_list = list(
#   obj1 = function(site_ids = c(), 
#                   pix_in_catch = c(), # a vector of all pixels across all catchments in design
#                   dist_mat = matrix(NA), 
#                   raster_stack = stack(), 
#                   pix_weights = c()){
#     # mean objective values of pix in catch (first stack in raster)
#     sum(raster_stack[[1]][pix_in_catch], na.rm=TRUE)
#   }
#   # option here for an objective involving power
# )
# 
# objective_list[[1]](pix_in_catch = site_catches[[1]],
#                     raster_stack = potential_vic_buffered)

# bad code !
# getNeighMat <- function(pixel_size, catchment_radius){
#   # assuming pixel_size and catchment_radius in same units ..
#   # produce neighbourhood matrix for call to adjacent()
#   mat_size = catchment_radius %/% pixel_size
#   
#   p <- st_as_sf(data.frame(x = 0,
#                            y = 0),
#                 coords=c("x","y"))
#   
#   out <- raster(matrix(1, nrow=mat_size*4 + 1, ncol=mat_size*4 + 1), 
#                 xmn=-(pixel_size*mat_size*2 + pixel_size/2), 
#                 xmx=pixel_size*mat_size*2 + pixel_size/2, 
#                 ymn=-(pixel_size*mat_size*2 + pixel_size/2), 
#                 ymx=pixel_size*mat_size*2 + pixel_size/2) %>% 
#     mask(st_buffer(p, catchment_radius)) %>%
#     trim() %>%
#     as.matrix()
#   
#   out[nrow(out)/2 + 1, ncol(out)/2 + 1] = 0
#   
#   out
# }


# neigh_mat = getNeighMat(pixel_size=res(potential_vic)[1],
#                         catchment_radius=0.05)


# nrow(site_ids) # lmao
# # This will take a while but only needs to be done once ..
# # like more than two hours BAD CODE !
# site_catches = lapply(site_ids[,"id"], function(x){adjacent(potential_vic_buffered, x, 
#                                                             directions=neigh_mat, 
#                                                             pairs=FALSE, 
#                                                             include=TRUE)})

# this needs to be some sort of product between the catchment membership matrix and the potential raster ..
# single_site_utility = sapply(1:nrow(site_ids), 
#                              function(x){objective_list[[1]](pix_in_catch = site_catches[[x]],
#                                                              raster_stack = potential_vic_buffered)})
single_site_utility <- matrix(potential_vic_buffered[catch_membership_mat],
                              nrow = nrow(catch_membership_mat),
                              ncol = ncol(catch_membership_mat)) %>%
  rowSums(na.rm = TRUE)

# check this matches what we expect :)
tmp <- potential_vic
tmp[!is.na(tmp)] <- single_site_utility

hist(single_site_utility, breaks=100)
# there are some zeroes in here ..
length(which(single_site_utility == 0)) # about 2%
single_site_quants = quantile(single_site_utility, probs=seq(0,1,0.05))

# visualise what we've done so far ...
plot(potential_vic_buffered, col=brewer.pal(9, "Reds"), xlim=c(140.9,141.2), ylim=c(-34.15,-33.9))
plot(potential_vic, col=brewer.pal(9,"Purples"), add=TRUE, new=TRUE, xlim=c(140.9,141.2), ylim=c(-34.15,-33.9))
plot(st_geometry(vic_shp), add=TRUE)


sel = single_site_utility > single_site_quants[20]
points(site_ids[sel, c("x","y")],
       col=alpha("orange", (single_site_utility[sel]-min(single_site_utility))/
                              (max(single_site_utility)-min(single_site_utility))),
       pch=16)


plot(potential_vic_buffered, col=brewer.pal(9, "Reds"))#, xlim=c(140.9,141.2), ylim=c(-34.15,-33.9))
plot(potential_vic, col=brewer.pal(9,"Purples"), add=TRUE, new=TRUE)#, xlim=c(140.9,141.2), ylim=c(-34.15,-33.9))
plot(st_geometry(vic_shp), add=TRUE)
points(site_ids[sel, c("x","y")],
       col=alpha("orange", (single_site_utility[sel]-min(single_site_utility))/
                   (max(single_site_utility)-min(single_site_utility))),
       pch=16, cex=0.1)

# rasterize single_site_utility?

# # are we ready to do some site selection? use distance-based proposal
# # suspect a broader interference radius will be necessary ..
# # this version can't even do 2*2 iters .. not sure where the complexity is coming from?
# # certainly need to set up for spartan ...
# # can i time individual steps and track cumulative energy?
# tmp = sim_an_anneal(site_catchment_list = site_catches,
#                     nselect = 10, 
#                     niters = 50, 
#                     raster_stack = potential_vic_buffered,
#                     obj_list = objective_list,
#                     find_union = TRUE,
#                     max_temperature = 1)
# # this is going to take a while ..
# # want heat map of frequently selected sites ..
# # still need to write some code to not look at designs we've already visited ..
# 
# # this is takes while for obvious reasons but I can't figure how to fix it :/
# ordered_single_site_utilities = order(single_site_utility, decreasing = TRUE)
# map_to_single_site_utility = sapply(1:length(single_site_utility), 
#                                     function(x){which(ordered_single_site_utilities == x)})
# # was meaning to chuck out that panel in favour of design utility anyway ...
# 
# wrap_sim_plots(tmp, panelled_plot = TRUE, single_site_utility = single_site_utility,
#                map_to_single_site_utility = 1:3,
#                main="Select 3 sites from sandbox: no annealing",
#                path = "~/Desktop/jev/output/site_select/start_jev_vic_mozzies.png")

# okey dokey ... I've got catch_membership_mat, let's fenangle it into sim_an_anneal
# require an objective_list

objective_list = list(
  obj1 = function(site_ids = c(), 
                  pix_in_catch = c(), # a vector of all pixels across all catchments in design - need to union() before calling objective
                  dist_mat = matrix(NA),
                  raster_stack = stack(),
                  pix_weights = c()){
    # mean objective values of pix in catch (first stack in raster)
    # if(is.na(sum(raster_stack[[1]][pix_in_catch]))){
    #   message("here")
    #   return(list(which(is.na(raster_stack[[1]][pix_in_catch]))))
    # }
    if (length(unique(raster_stack[[1]][pix_in_catch])) == 1){
      if(is.na(unique(raster_stack[[1]][pix_in_catch]))){message("here")}
  }
    sum(raster_stack[[1]][pix_in_catch], na.rm=TRUE)
  }
  # option here for an objective involving power
)

site_catches <- as.list(as.data.frame(t(catch_membership_mat)))
names(site_catches) <- site_ids[,"id"]

tmp <- potential_vic_buffered
tmp[is.na(tmp)] = 1
tmp <- as.data.frame(rasterToPoints(tmp))
tmp$id <- 1:nrow(tmp)
tmp2 <- tmp[which(tmp$id %in% unique(as.vector(catch_membership_mat))),]


compare_objectives = function(current, proposed, obj_list, site_catchment_list,
                              dist_mat, raster_stack, pix_weights, find_union){
  # given current and proposed sets of site IDs, report acceptance probability
  #message(current, proposed)
  # first, find union set of pixels in candidate surveillance design
  if (find_union == TRUE){
    current_net_catch = site_catchment_list[unlist(current)] %>%
      unlist() %>%
      unique()
    proposed_net_catch = site_catchment_list[unlist(proposed)] %>%
      unlist() %>%
      unique
  } else {
    current_net_catch = NA
    proposed_net_catch = NA
  }
  # for current/proposed, generate list of objective values
  current_objs = lapply(obj_list, function(f){f(current, current_net_catch,
                                                dist_mat, raster_stack,
                                                pix_weights)
  })
  
  #return(current_objs)
  proposed_objs = lapply(obj_list, function(f){f(proposed, proposed_net_catch,
                                                 dist_mat, raster_stack,
                                                 pix_weights)
  })
  # combine objective values
  pr_acc = apply(data.frame(unlist(current_objs),
                            unlist(proposed_objs)),
                 1, obj_to_pr) %>%
    prod()
  
  return(list(pr_acc = pr_acc,
              the_details = data.frame(unlist(current_objs),
                                       unlist(proposed_objs))))
}

compare_objectives(300, 170, objective_list, site_catches, 
                   raster_stack = raster::stack(potential_vic_buffered),
                   find_union = TRUE)


sim_an_anneal = function(site_catchment_list = list(), 
                         nselect = 0, 
                         niters = 0,
                         raster_stack = stack(),
                         dist_mat = matrix(NA),
                         site_catchment_weights_list = 0,  # hacking this for the poster
                         obj_list = list(),
                         find_union = FALSE,
                         max_temperature = 1,
                         init_sites = c()){
  
  nsites = length(site_catchment_list)
  # initialise object to keep track of selected sites
  outdf = data.frame(matrix(NA, 
                            nrow = niters*nselect + 1, 
                            ncol = nselect))
  if (length(init_sites) == nselect){
    outdf[1, ] = init_sites
  } else {
    outdf[1, ] = sample.int(nsites, 
                            size = nselect,
                            replace = FALSE)
  }
  
  temperatures = seq(1, max_temperature, length.out = nrow(outdf))
  acc_rat = 0
  pr_accs = c()
  deets = data.frame()
  step1 = 0
  step2 = 0
  step3 = 0
  # outer loop: niters (Metropolis)
  for (t in 1: niters){
    # inner loop: nselect (Gibbs)
    for (i in 1: nselect){
      # propose new design
      current = outdf[(t - 1) * nselect + i,]
      proposed = current
      t1 = Sys.time()
      # dropping this step as it's taking up a bunch of time?
      #probs = rep(1, nsites)
      #probs[unlist(current)] = 0
      proposed[i] = sample(x = as.vector(1: nsites)[-unlist(current)],
                           size = 1)#, # maybe prob is taking a while?
      #prob = sapply(1:nsites, function(x){ifelse(x %in% current, 0, 1)}))
      #prob = ifelse(1:nsites %in% current, 0, 1),
      #prob = probs)
      # message(paste0("current ", current, " proposed ", proposed))
      t2 = Sys.time()
      # calculate acceptance probability: obj(current_design, proposed_design)
      pr_acc = compare_objectives(current, proposed, obj_list,
                                  site_catchment_list,
                                  dist_mat, raster_stack, site_catchment_weights_list,
                                  find_union)
      # if(!"the_details" %in% names(pr_acc)){return(pr_acc)}
      t3 = Sys.time()
      pr_acc$pr_acc = pr_acc$pr_acc ** temperatures[(t - 1) * nselect + i]
      deets = rbind(deets, unlist(pr_acc$the_details))
      pr_acc = pr_acc$pr_acc
      pr_accs = c(pr_accs, unlist(pr_acc))
      t4 = Sys.time()
      
      step1 = step1 + t2 - t1
      step2 = step2 + t3 - t2
      step3 = step3 + t4 - t3
      # accept / reject
      if (runif(1) <= pr_acc){
        outdf[(t - 1) * nselect + i + 1,] = proposed
        acc_rat = acc_rat + 1
      } else {
        outdf[(t - 1) * nselect + i + 1,] = current
      }
    }
  }
  # names(deets) = as.vector(outer(c("Current", "Proposed"), 
  #                                1:length(obj_list), paste0))
  names(deets) = c(paste0("Current", 1:length(obj_list)),
                   paste0("Proposed", 1:length(obj_list)))
  
  return(list(outdf=outdf, # contains accepted designs
              acc_rat=acc_rat / (nrow(outdf) - 1),
              pr_accs=pr_accs,
              deets=deets, # contains utilities of current/proposed designs
              step1=step1,step2=step2,step3=step3))
}



tmp = sim_an_anneal(site_catchment_list = site_catches,
                    nselect = 10,
                    niters = 500,
                    raster_stack = raster::stack(potential_vic_buffered),
                    obj_list = objective_list,
                    find_union = TRUE,
                    max_temperature = 1)
# takes 30 seconds ...
30/500/60/60*1000000
# by my calculation, 16 hours for a million

# seems a weird spot ... need to check what's happening with buffering areas, site ids ...
# I'm too zonked to work out what's happening ...
# I have sussed the warnings and it's from the NaNs in the catch_membership_mat being used to index the raster
# so I'm okay with that

###############################################################################
# enumerate a little ...
{nrow(catch_membership_mat)
choose(nrow(catch_membership_mat), 10)
t1 = Sys.time()
# evaluate network distance objective also???
cand_designs <- sapply(1:10000000, function(x){sample(names(site_catches), 10, replace = FALSE)})
t2 = Sys.time()
# how many times were individual sites sampled?
unique(ftable(as.vector(cand_designs))[1,])
# how many sites weren't sampled at all?
length(which(!tmp %in% as.vector(cand_designs)))

length(unlist(site_catches[cand_designs[,2]]))
length(unique(unlist(site_catches[cand_designs[,2]])))

# takes a half hour or so ...
cand_shadows <- sapply(1:ncol(cand_designs), function(x){
  unlist(site_catches[cand_designs[,x]])
})
t3 = Sys.time()

dim(cand_shadows)

cand_shadows_unique <- sapply(1:ncol(cand_shadows), function(x){
  unique(cand_shadows[,x])
})
t4 = Sys.time()

design_utility <- lapply(cand_shadows_unique, function(x){
  sum(potential_vic_buffered[x], na.rm=TRUE)
})
t5 = Sys.time()

step1 = t2 - t1
step2 = step2 + t3 - t2
step3 = step3 + t4 - t3
step4 = step3 + t5 - t4
}

hist(unlist(design_utility), breaks=100)

design_utility_quants <- quantile(unlist(design_utility), probs=seq(0,1,0.00001))

abline(v=design_utility_quants, col="red")

tmp <- as.data.frame(rasterToPoints(potential_vic_buffered))
tmp$id <- which(!is.na(values(potential_vic_buffered)))

shortlist_designs <- cand_designs[,which(unlist(design_utility) > design_utility_quants[9990])]

tmp2 <- tmp[which(tmp$id %in% unlist(shortlist_designs)), c("x","y")]

plot(potential_vic_buffered, col=pn_cols)
points(tmp2, col=alpha("orange", 0.2), pch=16, cex=0.6)

shortlist_designs <- cand_designs[,which(unlist(design_utility) > design_utility_quants[100000])]
dim(shortlist_designs)
tmp2 <- tmp[which(tmp$id %in% unlist(shortlist_designs)), c("x","y")]

plot(potential_vic_buffered, col=pn_cols)
points(tmp2, col=alpha("orange", 0.4), pch=16, cex=0.6)

# still rather blotchy .. and that's the best 100 designs
