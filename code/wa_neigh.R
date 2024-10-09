setwd("/data/gpfs/projects/punim1228/jev_spartan/")
source("code/main.R")
source("code/genetic_algot.R")
source("code/wa_setup.R")

start_script <- Sys.time()

message("neigh 2")
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_neigh2 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_neigh2 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.2, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  png("figures/neigh2_test.png",
      height=1000,width=1000,pointsize=40)
  plot(raster(neigh_mat))
  plot(rasterToPolygons(raster(neigh_mat)), add=TRUE)
  dev.off()
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = wa_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = wa_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat,
                         pool = starting_point[[ind]], # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times_neigh2[ind] <- tend - tstart
    final_fronts_neigh2[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins for vic .. probably longer for WA .. but how much longer ?
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_neigh2_trapezoid.csv", row.names=FALSE)

save(times_neigh2,
     final_fronts_neigh2,
     file="output/diagnostics_wa_neigh_trapezoid.rds")

message(paste("neigh 3: ", Sys.time() - start_script))
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_neigh3 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_neigh3 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.3, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  png("figures/neigh3_test.png",
      height=1000,width=1000,pointsize=40)
  plot(raster(neigh_mat))
  plot(rasterToPolygons(raster(neigh_mat)), add=TRUE)
  dev.off()
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = wa_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = wa_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat,
                         pool = starting_point[[ind]], # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times_neigh3[ind] <- tend - tstart
    final_fronts_neigh3[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins for vic .. probably longer for WA .. but how much longer ?
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_neigh3_trapezoid.csv", row.names=FALSE)

save(times_neigh2, times_neigh3,
     final_fronts_neigh2, times_neigh3,
     file="output/diagnostics_wa_neigh_trapezoid.rds")

message(paste("neigh 4: ", Sys.time() - start_script))
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_neigh4 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_neigh4 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.4, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  png("figures/neigh4_test.png",
      height=1000,width=1000,pointsize=40)
  plot(raster(neigh_mat))
  plot(rasterToPolygons(raster(neigh_mat)), add=TRUE)
  dev.off()
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = wa_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = wa_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat,
                         pool = starting_point[[ind]], # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times_neigh4[ind] <- tend - tstart
    final_fronts_neigh4[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins for vic .. probably longer for WA .. but how much longer ?
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_neigh4_trapezoid.csv", row.names=FALSE)

save(times_neigh2, times_neigh3, times_neigh4,
     final_fronts_neigh2, final_fronts_neigh3, final_fronts_neigh4,
     file="output/diagnostics_wa_neigh_trapezoid.rds")

message(paste("neigh 5: ", Sys.time() - start_script))
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_neigh5 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_neigh5 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.45, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  png("figures/neigh5_test.png",
      height=1000,width=1000,pointsize=40)
  plot(raster(neigh_mat))
  plot(rasterToPolygons(raster(neigh_mat)), add=TRUE)
  dev.off()
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = wa_sites$id,
                         nselect = nselect, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = wa_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat,
                         pool = starting_point[[ind]], # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times_neigh5[ind] <- tend - tstart
    final_fronts_neigh5[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins for vic .. probably longer for WA .. but how much longer ?
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_neigh5_trapezoid.csv", row.names=FALSE)

save(times_neigh2, times_neigh3, times_neigh4, times_neigh5,
     final_fronts_neigh2, final_fronts_neigh3, final_fronts_neigh4, final_fronts_neigh5,
     file="output/diagnostics_wa_neigh_trapezoid.rds")


