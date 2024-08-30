# vic best settings?

setwd("/data/gpfs/projects/punim1228/jev_spartan/")
source("code/main.R")
source("code/genetic_algot.R")
source("code/vic_setup.R")

start_script <- Sys.time()

{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_apple <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_apple <- rep(list(NA), nruns)
  
  neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
  neigh_mat[!neigh_mat == 0] = 1
  catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
  catch_membership_mat <- values(catchment_stack, mat=TRUE)
  
  neigh_mat <- focalWeight(id_ras, 0.3, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 50000,
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
    message(paste("uninformed: ", ind, ";",tend-tstart))
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times_apple[ind] <- tend - tstart
    final_fronts_apple[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1}

write.csv(progress_auc, "output/vic_auc_pool50000_iters100_runs10_neigh3.csv", row.names=FALSE)
save(times_apple,
     final_fronts_apple,
     file="output/diagnostics_vic_50000_neigh3.rds")