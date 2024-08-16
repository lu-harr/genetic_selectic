# toy problem set up ....
getwd()

library(data.table)
library(rPref)

AGG_FACTOR = 10

guelphia_potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif') %>%
  crop(c(141, 146, -37.9, -32.9)) %>%
  aggregate(AGG_FACTOR) %>% # aggregated by mean .... just so you know ...
  t() %>%
  crop(c(-37.5, -36.65, 144.1, 144.9)) # don't ask :)

guelphia_hpop <- raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif') %>%
  crop(c(141, 146, -37.9, -32.9)) %>%
  aggregate(AGG_FACTOR) %>%
  t() %>%
  crop(c(-37.5, -36.65, 144.1, 144.9))

wtmat = matrix(1, nrow=3, ncol=3)
# should I worry about boundary effects? probably ...

toy_objective <- stack(guelphia_potential,
                             guelphia_hpop)
names(toy_objective) <- c("potent", "hpop")

id_ras <- toy_objective$potent
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
# now I want this as a list of vectors? Or a matrix?
catch_membership_mat <- values(catchment_stack, mat=TRUE)

site_ids <- data.frame(id=1:ncell(toy_objective),
                       potent=values(toy_objective$potent),
                       hpop=values(toy_objective$hpop))

# enumerate for all designs of five sites
# {t1 = Sys.time()
# enumerated <- combn(ncell(toy_objective$potent), 5, function(x){
#   c(x, 
#     sum(site_ids$hpop[unique(as.vector(catch_membership_mat[x,]))], na.rm=TRUE), 
#     sum(site_ids$potent[unique(as.vector(catch_membership_mat[x,]))], na.rm=TRUE))
# })
# t2 = Sys.time()}
# t2-t1 # 31 mins
# dim(enumerated)
# 
# {t1 = Sys.time()
# fwrite(t(enumerated),
#      file="~/Desktop/knowlesi/multi_site/output/toy_enumerated.csv")
# t2 = Sys.time()}
# it's 3 gigs!!

# {t1 = Sys.time()
enumerated <- fread("output/toy_enumerated.csv")
# t2 = Sys.time()}
names(enumerated) <- c(paste0("site", 1:5), "hpop", "potent")

# find exact pareto front
exact_toy_pareto <- psel(enumerated, high("potent")*high("hpop"))
exact_toy_pareto <- exact_toy_pareto[order(exact_toy_pareto$hpop)]

plot(exact_toy_pareto[,c("hpop","potent")])

pts_to_plot <- sample(nrow(enumerated), 200000)
eg_designs <- c(1,11,21)
toy_lonlats <- rasterToPoints(toy_objective)
eg_cols <- c("#c23375", "orange", "#8612ff")
eg_cex <- rev(seq(1.5,4,length.out=3))

get_catch_ras <- function(ras, ids){
  values(ras) <- NA
  values(ras)[site_ids$id[unique(as.vector(catch_membership_mat[ids,]))]] <- 1
  ras
}

# results figure: show exact solution
{png("figures/toy_problem_enumerated.png",
     height=1200, width=2000, pointsize=35)
  par(mfrow=c(1,2), mar=c(4.1,2.1,4.1,1.1), oma=c(0,2,0,2), xpd=NA)
  plot(enumerated[pts_to_plot, c("hpop","potent")],
       xlim=c(0, max(exact_toy_pareto$hpop)),
       ylim=c(1, max(exact_toy_pareto$potent)),
       xlab="Sum(Human Population Density)",
       ylab="Sum(Potential Risk)",
       main="Exact Pareto front")
  points(exact_toy_pareto[,c("hpop","potent")], col="black", pch=21, cex=1.3, bg=brat)
  points(exact_toy_pareto[eg_designs, c("hpop", "potent")],
         col="black", bg=eg_cols, pch=21, cex=1.6)
  
  # plot(guelphia_potential, col=pinks(100), axes=FALSE, type="n")
  par(mar=c(5.1,2.1,5.1,1.1))
  plot(rasterToPolygons(toy_objective$potent), border="black", col="white")

       #main="Examples of Pareto-optimal designs")
  #for (ind in 1:3){
  mtext("Examples of Pareto-optimal designs", side=3, line=2.5, cex=1.2, font=2)
  sites1 <- unlist(exact_toy_pareto[eg_designs[1],1:5])
  sites2 <- unlist(exact_toy_pareto[eg_designs[2],1:5])
  sites3 <- unlist(exact_toy_pareto[eg_designs[3],1:5])
  tmp1 <- get_catch_ras(toy_objective$potent, sites1)
  tmp2 <- get_catch_ras(toy_objective$potent, sites2)
  tmp3 <- get_catch_ras(toy_objective$potent, sites3)
  plot(tmp1, add=TRUE, col=alpha(eg_cols, 0.2)[1], legend=FALSE)
  plot(tmp2, add=TRUE, col=alpha(eg_cols, 0.2)[2], legend=FALSE)
  plot(tmp3, add=TRUE, col=alpha(eg_cols, 0.2)[3], legend=FALSE)
  # fiddling with cex to get them to overlap nicely
  points(toy_lonlats[sites1, c("x","y")], col=eg_cols[1], pch=15, cex=eg_cex[1])
  points(toy_lonlats[sites2, c("x","y")], col=eg_cols[2], pch=15, cex=eg_cex[c(1,2,1,2,1)])
  points(toy_lonlats[sites3, c("x","y")], col=eg_cols[3], pch=15, cex=eg_cex[c(2,3,1,2,1)])
  mtext("Potential risk \nhighest in north", 3, adj=0, col=eg_cols[1])
  mtext("Human population \nhighest in north-east", 3, adj=1, col=eg_cols[3])
  mtext("All catchments \noverlap", 1, adj=0.5, line=2)
  arrows(-37.43, 144.95, -37.35, 144.8, lwd=3, col=eg_cols[1])
  arrows(-36.7, 144.95, -36.75, 144.85, lwd=3, col=eg_cols[3])
  arrows(-37.05, 144.02, -37.01, 144.45, lwd=4)
  arrows(-37.05, 144.02, -36.85, 144.72, lwd=4)
    #lines(toy_lonlats[c(sites, sites[1]), c("x","y")], col=eg_cols[ind], lwd=2) 
    # an mst would be good here ...
    # don't quite care enough for now ...
  #}
  
  par(new=TRUE, oma=c(0,0,0,0), mar=c(0,0,0,0), mfrow=c(1,1))
  empty_plot_for_legend()
  subfigure_label(par()$usr, 0.1,0.93, "(a)", cex.label = 1.2)
  subfigure_label(par()$usr, 0.53,0.93, "(b)", cex.label = 1.2)
  
  dev.off()}

{png("figures/toy_problem_enumerated_pop_logged.png")
  par(mfrow=c(1,2))
  plot(log10(enumerated$hpop[pts_to_plot]),
       enumerated$potent[pts_to_plot],
       xlim=c(-0.5, log10(max(exact_toy_pareto$hpop))),
       ylim=c(1, max(exact_toy_pareto$potent)),
       xlab="log10(Human Population Density)",
       ylab="Potential Risk")
  points(log10(exact_toy_pareto$hpop),
         exact_toy_pareto$potent, col=brat, pch=16)
  dev.off()}


###############################################################################
# GA starts here

# let's see how we go with a random starting point ...

# also need to set box_extent

set.seed(834904)
starting_point <- matrix(sample(ncell(toy_objective), 2000, replace=TRUE), ncol=5) %>%
  as.data.frame()
names(starting_point) <- paste("site", 1:5, sep="")

educated_guess <- matrix(c(sample(order(values(toy_objective$potent), 
                                     decreasing = TRUE)[1:20], 1000, replace = TRUE),
                         sample(order(values(toy_objective$hpop),
                                      decreasing = TRUE)[1:20], 1000, replace = TRUE)),
                         ncol=5, byrow=TRUE) %>%
  as.data.frame()
names(educated_guess) <- paste("site", 1:5, sep="")

# let's run some experiments ...

{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  times50000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts50000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 50000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times50000[ind] <- tend - tstart
    final_fronts50000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool50000_iters100_runs10.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool50000_iters100_runs10.csv", row.names=FALSE)


{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  times10000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts10000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 10000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times10000[ind] <- tend - tstart
    final_fronts10000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # 2.5 mins per run ...
write.csv(progress_auc, "output/toy_auc_pool10000_iters100_runs10.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool10000_iters100_runs10.csv", row.names=FALSE)


{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times5000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts5000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 5000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times5000[ind] <- tend - tstart
    final_fronts5000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # 2.5 mins per run ...
write.csv(progress_auc, "output/toy_auc_pool5000_iters100_runs10.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool5000_iters100_runs10.csv", row.names=FALSE)


{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times1000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts1000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times1000[ind] <- tend - tstart
    final_fronts1000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
t2 = Sys.time()
t2-t1} #3.5 mins for 10 runs
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10.csv", row.names=FALSE)

# save(times1000, times5000, times10000, times50000,
#      final_fronts1000, final_fronts5000, final_fronts10000, final_fronts50000,
#      file="output/toy_diagnostics.rds")

progress_auc_1000 <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_auc_5000 <- read.csv("output/toy_auc_pool5000_iters100_runs10.csv")
progress_auc_10000 <- read.csv("output/toy_auc_pool10000_iters100_runs10.csv")
progress_auc_50000 <- read.csv("output/toy_auc_pool50000_iters100_runs10.csv")
# can restart by giving it current pool ...
# perhaps make it more explorative ?
# ... that would be change neighbourhood size .... from which we're sampling
# not sure how that scales?

matplot(progress_auc_10000[2:ncol(progress_auc_1000)], lty=1, col="black", type="l")
matplot(progress_auc_5000[2:ncol(progress_auc_1000)], lty=1, col="blue", type="l", add=TRUE)
matplot(progress_auc_1000[2:ncol(progress_auc_1000)], lty=1, col="red", type="l", add=TRUE)

auc_agg_fig(list(progress_auc_1000,
                  progress_auc_5000,
                  progress_auc_10000,
                  progress_auc_50000),
                 legend_labs=c("1,000", "5,000", "10,000", "50,000"),
                 legend_title="Pool size",
                 pal=iddu(4),
                 main="Increased pool size finds exact solution faster")

# would be nice if I could locate best run ... but let's just do a run
pareto_progress_contour(tmp$pareto_progress,
                       box_extent = c(95, max(exact_toy_pareto$hpop),
                                      3, max(exact_toy_pareto$potent)), # play around with this I guess ....
                       exact_soln = exact_toy_pareto[,c("hpop","potent")])


# how about an experiment with different numbers of neighbours included in each round?
# keep 1000 pool size from previous set of experiments
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  # much bigger than queen's case: (for sampling neighbours)
  neigh_mat <- focalWeight(id_ras, 0.2, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_stack <- mask(neigh_stack, rast(toy_objective$potent))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  neigh_membership_mat <- neigh_membership_mat[!is.na(values(toy_objective$potent)),]
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat, # second degree neighbours
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_neigh2.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_neigh2.csv", row.names=FALSE)

{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  # much bigger than queen's case: (for sampling neighbours)
  neigh_mat <- focalWeight(id_ras, 0.3, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_stack <- mask(neigh_stack, rast(toy_objective$potent))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  neigh_membership_mat <- neigh_membership_mat[!is.na(values(toy_objective$potent)),]
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat, # third degree neighbours
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_neigh3.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_neigh3.csv", row.names=FALSE)


progress_neigh1 <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_neigh2 <- read.csv("output/toy_auc_pool1000_iters100_runs10_neigh2.csv")
progress_neigh3 <- read.csv("output/toy_auc_pool1000_iters100_runs10_neigh3.csv")

auc_agg_fig(list(progress_neigh1[,2:ncol(progress_neigh1)],
                  progress_neigh2,
                  progress_neigh3),
             legend_labs=c("1st degree", "2nd degree", "3rd degree"),
             legend_title="Neighbourhood size",
             main="What do neighbours become?")
                           
                           
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  # much bigger than queen's case: (for sampling neighbours)
  neigh_mat <- focalWeight(id_ras, 0.2, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_stack <- mask(neigh_stack, rast(toy_objective$potent))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  neigh_membership_mat <- neigh_membership_mat[!is.na(values(toy_objective$potent)),]
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat, # second degree neighbours
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_neigh2.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_neigh2.csv", row.names=FALSE)


{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  # much bigger than queen's case: (for sampling neighbours)
  neigh_mat <- focalWeight(id_ras, 0.3, "circle")
  neigh_mat[!neigh_mat == 0] = 1
  neigh_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
  neigh_stack <- subset(neigh_stack, which(neigh_mat != 0))
  neigh_stack <- mask(neigh_stack, rast(toy_objective$potent))
  neigh_membership_mat <- values(neigh_stack, mat=TRUE)
  neigh_membership_mat <- neigh_membership_mat[!is.na(values(toy_objective$potent)),]
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = neigh_membership_mat, # third degree neighbours
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_neigh3.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_neigh3.csv", row.names=FALSE)

                           
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = educated_guess, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 1,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_loaded_start.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_loaded_start.csv", row.names=FALSE)

progress_uneducated <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_educated <- read.csv("output/toy_auc_pool1000_iters100_runs10_loaded_start.csv")

auc_agg_fig(list(progress_uneducated,
                  progress_educated),
                 legend_labs=c("Random", "Educated guess"),
                 legend_title="Starting pool",
                 main="Educated guess doesn't seem to help?")


# test for Pareto totality ?
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 2,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_pareto2.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_pareto2.csv", row.names=FALSE)

{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  progress_pc_pts <- matrix(NA, nrow=niters, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  
  for (ind in 1:nruns){
    tmp <- genetic_algot(site_ids = 1: nrow(site_ids),  # fix this - can be one function, but need to rewrite the same bit in the function
                         nselect = 5, 
                         poolsize = 1000,
                         niters = niters,
                         sandpit = toy_objective$potent,
                         potential_vec = site_ids$potent,
                         pop_vec = site_ids$hpop,
                         sample_method = "neighbours",
                         catchment_matrix = catch_membership_mat,
                         neighbourhood_matrix = catch_membership_mat, # keep it small for toy problem
                         pool = starting_point, # matrix of nselect columns
                         box_extent = c(0, max(exact_toy_pareto$hpop), 
                                        1, max(exact_toy_pareto$potent)),
                         top_level = 3,
                         plot_out = FALSE)
    
    progress_pc_pts[,ind] <- pareto_progress_pc_pts(tmp$pareto_progress, exact_toy_pareto)
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
  }
  
  t2 = Sys.time()
  t2-t1}
write.csv(progress_auc, "output/toy_auc_pool1000_iters100_runs10_pareto3.csv", row.names=FALSE)
write.csv(progress_pc_pts, "output/toy_pts_pool1000_iters100_runs10_pareto3.csv", row.names=FALSE)

progress_pareto1 <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_pareto2 <- read.csv("output/toy_auc_pool1000_iters100_runs10_pareto2.csv")
progress_pareto3 <- read.csv("output/toy_auc_pool1000_iters100_runs10_pareto3.csv")

auc_agg_fig(list(progress_pareto1,
                 progress_pareto2,
                 progress_pareto3),
            legend_labs=c("Rank 1", "Rank 1 & 2", "Rank 1, 2 & 3"),
            legend_title="Goldberg ranking",
            main="Making selection more lenient doesn't seem to change performance?")

################################################################################
# TURN THIS ALL INTO ONE FIGURE ... TWEAK COLOURS ...

mean(times1000)
sd(times1000)
mean(times5000)
sd(times5000)
mean(times10000)
sd(times10000)
mean(times50000)
sd(times50000)

alltimes <- data.frame(poolsize=c(rep(1000, 10), rep(5000, 10), rep(10000, 10), rep(50000, 10)),
                       time=c(times1000, times5000*60, times10000*60, times50000*60))
timemod <- lm(time ~ poolsize, alltimes)
summary(timemod)
plot(alltimes)


progress_auc_1000 <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_auc_5000 <- read.csv("output/toy_auc_pool5000_iters100_runs10.csv")
progress_auc_10000 <- read.csv("output/toy_auc_pool10000_iters100_runs10.csv")
progress_auc_50000 <- read.csv("output/toy_auc_pool50000_iters100_runs10.csv")

progress_neigh1 <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_neigh2 <- read.csv("output/toy_auc_pool1000_iters100_runs10_neigh2.csv")
progress_neigh3 <- read.csv("output/toy_auc_pool1000_iters100_runs10_neigh3.csv")

progress_uneducated <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_educated <- read.csv("output/toy_auc_pool1000_iters100_runs10_loaded_start.csv")

progress_pareto1 <- read.csv("output/toy_auc_pool1000_iters100_runs10.csv")
progress_pareto2 <- read.csv("output/toy_auc_pool1000_iters100_runs10_pareto2.csv")
progress_pareto3 <- read.csv("output/toy_auc_pool1000_iters100_runs10_pareto3.csv")

# Oh no! Now I need to un-subtract exact AUC from all the progress_paretooooooss!

exact_soln <- as.data.frame(exact_toy_pareto)[,grep("site", names(exact_toy_pareto), invert=TRUE)]
exact_soln <- exact_soln[order(exact_soln[,1]),]
exact_soln <- rbind(c(0, max(exact_soln[,2])),
                    exact_soln,
                    c(max(exact_soln[,1]), 0))
exact_auc <- area_under_curve(exact_soln[,1], exact_soln[,2], method="step")

{png("figures/toy_sensitivity.png",
     width=2400,
     height=2000,
     pointsize = 40)

par(mfrow=c(2,2), mar=c(2.1,2.1,2.1,2.1), oma=c(3,4,1,0))

auc_agg_fig(list(progress_auc_1000,
                progress_auc_5000,
                progress_auc_10000,
                progress_auc_50000),
             legend_labs=c("1,000", "5,000", "10,000", "50,000"),
             legend_title="Pool size",
              ylim = c(890, exact_auc),
              pal=c(iddu(4)[2:4], brat))#,
                 #main="Increased pool size finds exact solution faster")

auc_agg_fig(list(progress_neigh1,
                 progress_neigh2,
                 progress_neigh3),
            legend_labs=c("1st degree", "2nd degree", "3rd degree"),
            legend_title="Neighbourhood size",
            ylim = c(890, exact_auc),
            pal=iddu(4)[2:4])
                  #main="What do neighbours become?")

auc_agg_fig(list(progress_pareto1,
                 progress_pareto2,
                 progress_pareto3),
            legend_labs=c("Rank 1", "Rank 1 & 2", "Rank 1, 2 & 3"),
            legend_title="Goldberg ranking",
            ylim = c(890, exact_auc),
            pal=iddu(4)[2:4])
            #main="Making selection more lenient doesn't seem to change performance?")

auc_agg_fig(list(progress_uneducated,
                 progress_educated),
            legend_labs=c("Random", "Educated guess"),
            legend_title="Starting pool",
            ylim = c(890, exact_auc),
            pal=iddu(4)[2:4])
#main="Educated guess doesn't seem to help?")

mtext("Area under estimated Pareto front", 2, outer=TRUE, line=2, cex=1.2)
mtext("Iteration", 1, outer=TRUE, line=1, cex=1.2)

par(new=TRUE, mfrow=c(1,1), xpd=NA)
empty_plot_for_legend()
subfigure_label(par()$usr, -0.05, 1.06, "(a)")
subfigure_label(par()$usr, 0.515, 1.06, "(b)")
subfigure_label(par()$usr, -0.05, 0.48, "(c)")
subfigure_label(par()$usr, 0.515, 0.48, "(d)")
dev.off()}

# looks like I got a bunch of those points !
final_frontsdf <- rbindlist(final_fronts10000) %>%
  as.data.frame()
tmp <- final_frontsdf[,grep("site", names(final_frontsdf))]
tmp <- as.data.frame(t(apply(tmp, 1, sort)))
names(tmp) <- paste("site", 1:ncol(tmp), sep="")
final_frontsdf[,grep("site", names(final_frontsdf))] <- tmp
collected <- final_frontsdf %>%
  dplyr::select(grep("site", names(final_frontsdf), value=TRUE)) %>%
  # there appears to be some mismatch ... is it a problem with the ordering step ??
  # might be a rounding problem in the objective columns? disappears when I unique over site columns ...
  unique() %>% # only 100 now
  inner_join(exact_toy_pareto, by = grep("site", names(final_frontsdf), value=TRUE))

# contour plot but now it's final fronts only .. check this for 10,000 iters
# work out a way to save everything? might just have to be an rds
pal = viridis(length(final_fronts10000) + 2)[1:length(final_fronts10000)]
pal=rep(viridis(1), length(final_fronts10000))
#plot(0, xlim=c(0, max(exact_toy_pareto$hpop)), ylim=c(1, max(exact_toy_pareto$potent)), type="n")
plot(0, xlim=c(100, max(exact_toy_pareto$hpop)), 
     ylim=c(4.5, max(exact_toy_pareto$potent)), type="n",
     main="14/21 points recovered in 25 minutes computing time",
     xlab="Sum(Human pop)", ylab="Sum(Potential risk)")
for (ind in 1:length(final_fronts10000)){
  front = final_fronts10000[[ind]][,c("sum_pop", "sum_risk")]
  front <- front[order(front$sum_pop),]
  points(front, col=pal[ind])
  lines(front, col=pal[ind])
}
points(exact_toy_pareto[,c("hpop", "potent")], col=iddu(4)[4], pch=16)
lines(exact_toy_pareto[,c("hpop", "potent")], col=iddu(4)[4])
points(collected[,c("hpop", "potent")], pch=16, col=brat)
# niceee
