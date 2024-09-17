setwd("/data/gpfs/projects/punim1228/jev_spartan/")
source("code/main.R")
source("code/genetic_algot.R")
source("code/vic_setup.R")
# note that this takes three goes as is :)

start_script <- Sys.time()

# # pool size 1000
# message("pool size 1,000")
# 
# {set.seed(834903)
#   t1 = Sys.time()
#   niters = 100
#   nruns = 10
# 
#   times1000 <- matrix(NA, nrow=1, ncol=nruns)
#   progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
#   final_fronts1000 <- rep(list(NA), nruns)
# 
#   neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
#   neigh_mat[!neigh_mat == 0] = 1
#   catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
#   catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
#   catch_membership_mat <- values(catchment_stack, mat=TRUE)
# 
#   for (ind in 1:nruns){
#     tstart <- Sys.time()
#     tmp <- genetic_algot(site_ids = vic_sites$id,
#                          nselect = nselect,
#                          poolsize = 1000,
#                          niters = niters,
#                          sandpit = vic_objective$potent,
#                          potential_vec = site_ids$potent,
#                          pop_vec = site_ids$hpop,
#                          sample_method = "neighbours",
#                          catchment_matrix = catch_membership_mat,
#                          neighbourhood_matrix = catch_membership_mat,
#                          pool = starting_point, # matrix of nselect columns
#                          top_level = 1,
#                          plot_out = FALSE)
#     tend <- Sys.time()
# 
#     progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
#     times1000[ind] <- tend - tstart
#     final_fronts1000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
#   }
# 
#   t2 = Sys.time()
#   t2-t1}
# 
# write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10_trapezoid.csv", row.names=FALSE)
# save(times1000,
#      final_fronts1000,
#      file="output/diagnostics_vic_pool_trapezoid.rds")
# 
# 
# message(paste("pool size 5,000: ", Sys.time() - start_script))
# {set.seed(834903)
#   t1 = Sys.time()
#   niters = 100
#   nruns = 10
# 
#   times5000 <- matrix(NA, nrow=1, ncol=nruns)
#   progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
#   final_fronts5000 <- rep(list(NA), nruns)
# 
#   neigh_mat <- focalWeight(id_ras, 0.12, "circle")
#   neigh_mat[!neigh_mat == 0] = 1
#   catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
#   catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
#   catch_membership_mat <- values(catchment_stack, mat=TRUE)
# 
#   for (ind in 1:nruns){
#     tstart <- Sys.time()
#     tmp <- genetic_algot(site_ids = vic_sites$id,
#                          nselect = nselect,
#                          poolsize = 5000,
#                          niters = niters,
#                          sandpit = vic_objective$potent,
#                          potential_vec = site_ids$potent,
#                          pop_vec = site_ids$hpop,
#                          sample_method = "neighbours",
#                          catchment_matrix = catch_membership_mat,
#                          neighbourhood_matrix = catch_membership_mat,
#                          pool = starting_point,
#                          top_level = 1,
#                          plot_out = FALSE)
#     tend <- Sys.time()
#     message(paste("pool size 5,000: ", ind, ";",tend-tstart))
#     progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
#     times5000[ind] <- tend - tstart
#     final_fronts5000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
#   }
# 
#   t2 = Sys.time()
#   t2-t1}
# 
# write.csv(progress_auc, "output/vic_auc_pool5000_iters100_runs10_trapezoid.csv", row.names=FALSE)
# save(times1000, times5000,
#      final_fronts1000, final_fronts5000,
#      file="output/diagnostics_vic_pool_trapezoid.rds")
# 
# 
# message(paste("pool size 10,000: ", Sys.time() - start_script))
# {set.seed(834903)
#   t1 = Sys.time()
#   niters = 100
#   nruns = 10
# 
#   times10000 <- matrix(NA, nrow=1, ncol=nruns)
#   progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
#   final_fronts10000 <- rep(list(NA), nruns)
# 
#   neigh_mat <- focalWeight(id_ras, 0.12, "circle")
#   neigh_mat[!neigh_mat == 0] = 1
#   catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
#   catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
#   catch_membership_mat <- values(catchment_stack, mat=TRUE)
# 
#   for (ind in 1:nruns){
#     tstart <- Sys.time()
#     tmp <- genetic_algot(site_ids = vic_sites$id,
#                          nselect = nselect,
#                          poolsize = 10000,
#                          niters = niters,
#                          sandpit = vic_objective$potent,
#                          potential_vec = site_ids$potent,
#                          pop_vec = site_ids$hpop,
#                          sample_method = "neighbours",
#                          catchment_matrix = catch_membership_mat,
#                          neighbourhood_matrix = catch_membership_mat,
#                          pool = starting_point,
#                          top_level = 1,
#                          plot_out = FALSE)
#     tend <- Sys.time()
#     message(paste("pool size 10,000: ", ind, ";",tend-tstart))
#     progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
#     times10000[ind] <- tend - tstart
#     final_fronts10000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
#   }
# 
#   t2 = Sys.time()
#   t2-t1} # expecting 15*5 mins
# write.csv(progress_auc, "output/vic_auc_pool10000_iters100_runs10_trapezoid.csv", row.names=FALSE)
# save(times1000, times5000, times10000,
#      final_fronts1000, final_fronts5000, final_fronts10000,
#      file="output/diagnostics_vic_pool_trapezoid.rds")
# 
# load(file="output/diagnostics_vic_pool_trapezoid.rds")
# 
# message(paste("pool size 50,000: ", Sys.time() - start_script))
# {set.seed(834903)
#   t1 = Sys.time()
#   niters = 100
#   nruns = 10
# 
#   times50000 <- matrix(NA, nrow=1, ncol=nruns)
#   progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
#   final_fronts50000 <- rep(list(NA), nruns)
# 
#   neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
#   neigh_mat[!neigh_mat == 0] = 1
#   catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
#   catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
#   catch_membership_mat <- values(catchment_stack, mat=TRUE)
# 
#   for (ind in 1:nruns){
#     tstart <- Sys.time()
#     tmp <- genetic_algot(site_ids = vic_sites$id,
#                          nselect = nselect,
#                          poolsize = 50000,
#                          niters = niters,
#                          sandpit = vic_objective$potent,
#                          potential_vec = site_ids$potent,
#                          pop_vec = site_ids$hpop,
#                          sample_method = "neighbours",
#                          catchment_matrix = catch_membership_mat,
#                          neighbourhood_matrix = catch_membership_mat,
#                          pool = starting_point, # matrix of nselect columns
#                          top_level = 1,
#                          plot_out = FALSE)
#     tend <- Sys.time()
#     message(paste("pool size 50,000: ", ind, ";",tend-tstart))
#     progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
#     times50000[ind] <- tend - tstart
#     final_fronts50000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
#   }
# 
#   t2 = Sys.time()
#   t2-t1}
# write.csv(progress_auc, "output/vic_auc_pool50000_iters100_runs10_trapezoid.csv", row.names=FALSE)
# 
# save(times1000, times5000, times10000, times50000,
#      final_fronts1000, final_fronts5000, final_fronts10000, final_fronts50000,
#      file="output/diagnostics_vic_pool_trapezoid.rds")

load(file="output/diagnostics_vic_pool_trapezoid.rds")

message(paste("pool size 100,000: ", Sys.time() - start_script))
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times100000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts100000 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
  neigh_mat[!neigh_mat == 0] = 1
  catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
  catch_membership_mat <- values(catchment_stack, mat=TRUE)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 100000,
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
    message(paste("pool size 100,000: ", ind, ";",tend-tstart))
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times50000[ind] <- tend - tstart
    final_fronts50000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/vic_auc_pool100000_iters100_runs10_trapezoid.csv", row.names=FALSE)

save(times1000, times5000, times10000, times50000, times100000,
     final_fronts1000, final_fronts5000, final_fronts10000, final_fronts50000, final_fronts100000,
     file="output/diagnostics_vic_pool_trapezoid.rds")