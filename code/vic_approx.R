vic_shp <- states %>%
  filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)

# need 20km to do third-degree neighbours ...
vic_shadow <- vic_shp %>%
  st_buffer(dist = 25000) %>% # 80% certain we're in metres
  st_simplify(dTolerance = 1000)
# ugly ! ah well !

vic_objective <- stack(raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif'),
                       raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif')) %>%
                         aggregate(AGG_FACTOR) %>%
                         crop(extent(vic_shadow)) %>%
                         mask(vic_shadow)
names(vic_objective) <- c("potent", "hpop")

tmp <- vic_objective
vic_objective <- mask(vic_objective, vic_shp)
vic_objective$buffer_potent <- tmp$potent
vic_objective$buffer_hpop <- tmp$hpop

id_ras <- vic_objective$potent
id_ras[] <- 1:ncell(id_ras)

# I guess compare performances by plotting fronts next to each other?
# Stop when there hasn't been a change for 10 iters?

# try unaggregated case first ...

# raster level assessment of existing trapping!
# all_mozzies is read in in sa4_appraisal ... should transport into main.R

# just want unique locations
vic_mozzies <- subset(all_mozzies, data_source == "Vic") %>%
  filter(as.Date(date) > as.Date("2022-07-01")) %>% 
  # consider only 2022-2023 trapping season
  group_by(longitude, latitude) %>%
  summarise(ntimes = n()) %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise()
# now need to re-assign latlons ... 

# hist(as.Date(vic_mozzies$date), 
#      breaks=as.Date(c(paste0("2022-", 7:12, "-01"), paste0("2023-", 1:7, "-01"))), 
#      freq=TRUE)
  
# raster::extract(id_ras, vic_mozzies[,c("longitude", "latitude")])
  


neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
# don't actually want *weights* per se .. but that might be useful later ...
neigh_mat[!neigh_mat == 0] = 1
# jumping into terra because of option for fun argument that returns vector
# this gives me a layered SpatRaster
catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c) # so much quicker than my diy version
# need to remove zero layers ...
catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
# and remove non-victorian pixels that were part of the buffer ...
# now I want this as a list of vectors? Or a matrix?
catch_membership_mat <- values(catchment_stack, mat=TRUE)

# making an edit here to use buffered values? As these end up in catchment calculations?
site_ids <- as.data.frame(rasterToPoints(id_ras))
site_ids$potent <- values(vic_objective$buffer_potent)
site_ids$hpop <- values(vic_objective$buffer_hpop)
site_ids$id <- 1:ncell(vic_objective)

vic_sites <- site_ids[!is.na(values(vic_objective$potent)),]
# 3297 victorian sites to choose from ... which I guess halves the problem


################################################################################
# let's deploy GA and see what happens
nselect <- 100

set.seed(834904)
starting_point <- matrix(sample(vic_sites$id, 100000, replace=TRUE), ncol=nselect) %>%
  as.data.frame()
names(starting_point) <- paste("site", 1:nselect, sep="")

educated_guess <- matrix(c(sample(vic_sites$id[order(vic_sites$potent, 
                                                     decreasing = TRUE)][1:200], 50000, replace = TRUE),
                           sample(vic_sites$id[order(vic_sites$hpop, 
                                               decreasing = TRUE)][1:200], 50000, replace = TRUE)),
                         ncol=nselect, byrow=TRUE) %>%
  as.data.frame()
names(educated_guess) <- paste("site", 1:nselect, sep="")

# would be nice if I could plot the educated guess to assess coverage ...
plot(vic_objective$potent)
for (i in 1:nrow(educated_guess)){
  tmp <- which(vic_sites$id %in% educated_guess[i,])
  points(vic_sites[tmp, c("x","y")])
}
# that takes a sec but is what I expected ... 
# perhaps try starting from a risk-only point later?

niters = 100
{set.seed(834903)
  tstart1 <- Sys.time()
tmp <- genetic_algot(site_ids = vic_sites$id,
                     nselect = nselect, 
                     poolsize = 1000,
                     niters = niters,
                     sandpit = vic_objective$potent,
                     potential_vec = site_ids$potent,
                     pop_vec = site_ids$hpop,
                     sample_method = "neighbours",
                     catchment_matrix = catch_membership_mat,
                     neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                     pool = starting_point, # matrix of nselect columns
                     top_level = 1,
                     shp = vic_shp,
                     box_extent = c(5000,45000,35,90),
                     plot_out = TRUE)
tend1 <- Sys.time()} # 1.7 mins to do 100 iters

# baseline run: poolsize 1000

{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times1000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts1000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = vic_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat,
                         pool = starting_point, # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times1000[ind] <- tend - tstart
    final_fronts1000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins
write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10.csv", row.names=FALSE)


{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times5000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts5000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 5000,
                         niters = niters,
                         sandpit = vic_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat,
                         pool = starting_point, # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times5000[ind] <- tend - tstart
    final_fronts5000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15*5 mins
write.csv(progress_auc, "output/vic_auc_pool5000_iters100_runs10.csv", row.names=FALSE)


auc_agg_fig(list(progress_auc)) # interesting

################################################################################
# varying neighbourhood size experiment ... 1,2,3 might be too small ... see if times differ and try 10 at some point

# including second-degree neighbours
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times1000neigh2 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts1000neigh2 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.2, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = vic_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat,
                         pool = starting_point, # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times1000neigh2[ind] <- tend - tstart
    final_fronts1000neigh2[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins
write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10_neigh2.csv", row.names=FALSE)

# including third-degree neighbours
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times1000neigh3 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts1000neigh3 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.3, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  for (ind in 1:nruns){
    message(ind)
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = vic_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat,
                         pool = starting_point, # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times1000neigh3[ind] <- tend - tstart
    final_fronts1000neigh3[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins
write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10_neigh3.csv", row.names=FALSE)



progress_neigh1 <- read.csv("output/vic_auc_pool1000_iters100_runs10.csv")
progress_neigh2 <- read.csv("output/vic_auc_pool1000_iters100_runs10_neigh2.csv")
progress_neigh3 <- read.csv("output/vic_auc_pool1000_iters100_runs10_neigh3.csv")

# very interesting ...
auc_agg_fig(list(progress_neigh1,
                 progress_neigh2),
            legend_labs=c("1st degree", "2nd degree", "3rd degree"),
            legend_title="Neighbourhood size",
            main="What do neighbours become?")
times1000neigh2


  
  






