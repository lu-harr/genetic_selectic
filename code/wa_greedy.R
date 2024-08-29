setwd("/data/gpfs/projects/punim1228/jev_spartan/")
source("code/main.R")
source("code/genetic_algot.R")
source("code/wa_setup.R")

start_script <- Sys.time()

message("here we go: greedy start")
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_greedy <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_greedy <- rep(list(NA), nruns)
  
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
                         pool = educated_guess, # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times_greedy[ind] <- tend - tstart
    final_fronts_greedy[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins for vic .. probably longer for WA .. but how much longer ?
write.csv(progress_auc, "output/wa_auc_pool1000_iters100_runs10_greedystart.csv", row.names=FALSE)

save(times_greedy,
     final_fronts_greedy,
     file="output/diagnostics_wa_greedy.rds")

