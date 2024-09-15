setwd("/data/gpfs/projects/punim1228/jev_spartan/")
source("code/main.R")
source("code/genetic_algot.R")
source("code/wa_setup.R")

start_script <- Sys.time()

message("pareto 2")
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_pareto2 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_pareto2 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
  neigh_mat[!neigh_mat == 0] = 1
  catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
  catch_membership_mat <- values(catchment_stack, mat=TRUE)
  
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
                         neighbourhood_matrix = catch_membership_mat,
                         pool = starting_point, # matrix of nselect columns
                         top_level = 2,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times_pareto2[ind] <- tend - tstart
    final_fronts_pareto2[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_pareto2_trapezoid.csv", row.names=FALSE)

save(times_pareto2,
     final_fronts_pareto2,
     file="output/diagnostics_wa_pareto_trapezoid.rds")

message(paste("pareto 3: ", Sys.time() - start_script))
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_pareto3 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_pareto3 <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
  neigh_mat[!neigh_mat == 0] = 1
  catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
  catch_membership_mat <- values(catchment_stack, mat=TRUE)
  
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
                         neighbourhood_matrix = catch_membership_mat,
                         pool = starting_point, # matrix of nselect columns
                         top_level = 3,
                         plot_out = FALSE)
    tend <- Sys.time()
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress, method="trapezoid")
    times_pareto3[ind] <- tend - tstart
    final_fronts_pareto3[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_pareto3_trapezoid.csv", row.names=FALSE)

save(times_pareto2, times_pareto3,
     final_fronts_pareto2, final_fronts_pareto3,
     file="output/diagnostics_wa_pareto_trapezoid.rds")


