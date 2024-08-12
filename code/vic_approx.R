# all of the analysis and figures for Victorian problem are in here
# (calls to genetic_algot.R)

vic_shp <- states %>%
  filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)

# need 20km to do third-degree neighbours ...
vic_shadow <- vic_shp %>%
  st_buffer(dist = 25000) %>% # 80% certain we're in metres
  st_simplify(dTolerance = 1000)
# ugly ! ah well !

AGG_FACTOR <- 10
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
# 124 left

# hist(as.Date(vic_mozzies$date), 
#      breaks=as.Date(c(paste0("2022-", 7:12, "-01"), paste0("2023-", 1:7, "-01"))), 
#      freq=TRUE)

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


plot(vic_objective$potent)
#points(vic_mozzies$longitude, vic_mozzies$latitude, pch=16)
points(site_ids[vic_mozzies$pix, c("x","y")], col="blue", pch=16)
objective_func = function(x, catch_mem, vec){
  sum(vec[unique(as.vector(catch_mem[unlist(x),]))], na.rm=TRUE)
}

# not sure if I believe how good this is under risk obj ... although there are more sites here
existing_potent <- objective_func(vic_mozzies$pix, 
                                  catch_membership_mat, 
                                  vic_objective$buffer_potent)
existing_hpop <- objective_func(vic_mozzies$pix,
                                catch_membership_mat,
                                vic_objective$buffer_hpop)

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

################################################################################
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


{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times10000 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts10000 <- rep(list(NA), nruns)
  
  for (ind in 1:nruns){
    tstart <- Sys.time()
    tmp <- genetic_algot(site_ids = vic_sites$id,
                         nselect = nselect, 
                         poolsize = 10000,
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
    times10000[ind] <- tend - tstart
    final_fronts10000[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15*5 mins
write.csv(progress_auc, "output/vic_auc_pool10000_iters100_runs10.csv", row.names=FALSE)


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

################################################################################
# educated guess
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_educated <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_educated <- rep(list(NA), nruns)
  
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
                         pool = educated_guess, # matrix of nselect columns
                         top_level = 1,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times_educated[ind] <- tend - tstart
    final_fronts_educated[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins
write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10_loaded_start.csv", row.names=FALSE)

# want to see the turn-around?
# would be good if I could report AUF .... just deploy contour plot
tmp <- genetic_algot(site_ids = vic_sites$id,
                     nselect = nselect, 
                     poolsize = 1000,
                     niters = 15,
                     sandpit = vic_objective$potent,
                     potential_vec = site_ids$potent,
                     pop_vec = site_ids$hpop,
                     sample_method = "neighbours",
                     catchment_matrix = catch_membership_mat,
                     neighbourhood_matrix = catch_membership_mat,
                     pool = educated_guess, # matrix of nselect columns
                     top_level = 1,
                     plot_out = TRUE)

tmp3 = tmp$pareto_progress[[1]][,c("sum_pop", "sum_risk")]
tmp3 <- rbind(c(0, max(tmp3$sum_risk)),
              tmp3[order(tmp3$sum_pop),],
              c(max(tmp3$sum_pop), 0))
tmp4 = tmp$pareto_progress[[2]][,c("sum_pop", "sum_risk")]
tmp4 <- rbind(c(0, max(tmp4$sum_risk)),
              tmp4[order(tmp4$sum_pop),],
              c(max(tmp4$sum_pop), 0))

tmp2 = pareto_progress_contour(list(tmp3[,2:1], tmp4[,2:1]), plot_auc=TRUE)#, 
                              # box_extent=c(15000, 45000, 30, 105))



area_under_curve(tmp3$sum_risk[1:(nrow(tmp3)-1)],
                 tmp3$sum_pop[1:(nrow(tmp3)-1)],
                 method="step")
area_under_curve(tmp$pareto_progress[[1]][,"sum_risk"],
                 tmp$pareto_progress[[1]][,"sum_pop"],
                 method="step")
area_under_curve(tmp$pareto_progress[[2]][,"sum_risk"],
                 tmp$pareto_progress[[2]][,"sum_pop"],
                 method="step")
# definitely do need to do that rbind ...

##################################################################################
# Pareto optimality
# 11:47
{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_pareto2 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_pareto2 <- rep(list(NA), nruns)
  
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
                         top_level = 2,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times_pareto2[ind] <- tend - tstart
    final_fronts_pareto2[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins
write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10_pareto2.csv", row.names=FALSE)

{set.seed(834903)
  t1 = Sys.time()
  niters = 100
  nruns = 10
  
  times_pareto3 <- matrix(NA, nrow=1, ncol=nruns)
  progress_auc <- matrix(NA, nrow=niters, ncol=nruns)
  final_fronts_pareto3 <- rep(list(NA), nruns)
  
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
                         top_level = 3,
                         plot_out = FALSE)
    tend <- Sys.time()
    
    progress_auc[,ind] <- pareto_progress_auc(tmp$pareto_progress)
    times_pareto3[ind] <- tend - tstart
    final_fronts_pareto3[[ind]] <- tmp$pareto_progress[[length(tmp$pareto_progress)]]
  }
  
  t2 = Sys.time()
  t2-t1} # expecting 15 mins
write.csv(progress_auc, "output/vic_auc_pool1000_iters100_runs10_pareto3.csv", row.names=FALSE)


##################################################################################
# LET'S REVIEW ...

progress_neigh1 <- read.csv("output/vic_auc_pool1000_iters100_runs10.csv")
progress_neigh2 <- read.csv("output/vic_auc_pool1000_iters100_runs10_neigh2.csv")
progress_neigh3 <- read.csv("output/vic_auc_pool1000_iters100_runs10_neigh3.csv")

progress1000 <- read.csv("output/vic_auc_pool1000_iters100_runs10.csv")
progress5000 <- read.csv("output/vic_auc_pool5000_iters100_runs10.csv")
progress10000 <- read.csv("output/vic_auc_pool10000_iters100_runs10.csv")

progress_uneducated <- read.csv("output/vic_auc_pool1000_iters100_runs10.csv")
progress_educated <- read.csv("output/vic_auc_pool1000_iters100_runs10_loaded_start.csv")

progress_pareto1 <- read.csv("output/vic_auc_pool1000_iters100_runs10.csv")
progress_pareto2 <- read.csv("output/vic_auc_pool1000_iters100_runs10_pareto2.csv")
progress_pareto3 <- read.csv("output/vic_auc_pool1000_iters100_runs10_pareto3.csv")

save(times1000, times5000, times10000,
     times1000neigh2, times1000neigh3,
     times_educated,
     times_pareto2, times_pareto3,
     final_fronts1000, final_fronts5000, final_fronts10000,
     final_fronts1000neigh2, final_fronts1000neigh3,
     final_fronts_educated,
     final_fronts_pareto2,final_fronts_pareto3,
     file="output/vic_diagnostics.rds")

load("output/vic_diagnostics.rds")

# very interesting ...
# might need to make the box_extents the same for all of these for context ...?

{png("figures/vic_sensitivity.png",
     width=2400,
     height=2000,
     pointsize = 40)
  
  par(mfrow=c(2,2), mar=c(2.1,2.1,2.1,2.1), oma=c(3,4,1,0))

auc_agg_fig(list(progress1000,
                 progress5000,
                 progress10000),
            legend_labs=c("1,000", "5,000", "10,000"),
            legend_title="Pool size",
            pal=iddu(4)[2:4])

# a bit concerned that this line goes down ....... investigate ....
# could be because of concavitity but doubt it happens as much as appears?

auc_agg_fig(list(progress_neigh1,
                 progress_neigh2,
                 progress_neigh3),
            legend_labs=c("1st degree", "2nd degree", "3rd degree"),
            legend_title="Neighbourhood size",
            pal=iddu(4)[2:4])

auc_agg_fig(list(progress_pareto1,
                 progress_pareto2,
                 progress_pareto3),
            legend_labs=c("Rank 1", "Rank 1 & 2", "Rank 1, 2 & 3"),
            legend_title="Goldberg ranking",
            pal=iddu(4)[2:4])

auc_agg_fig(list(progress_uneducated,
                 progress_educated),
            legend_labs=c("Random", "Educated guess"),
            legend_title="Starting pool",
            pal=iddu(4)[2:4])

mtext("Area under estimated Pareto front", 2, outer=TRUE, line=2, cex=1.2)
mtext("Iteration", 1, outer=TRUE, line=1, cex=1.2)

par(new=TRUE, mfrow=c(1,1), xpd=NA)
empty_plot_for_legend()
subfigure_label(par()$usr, -0.05, 1.06, "(a)")
subfigure_label(par()$usr, 0.515, 1.06, "(b)")
subfigure_label(par()$usr, -0.05, 0.48, "(c)")
subfigure_label(par()$usr, 0.515, 0.48, "(d)")
dev.off()}

mean(times1000)
mean(times1000neigh2)
mean(times1000neigh3) # does take a tiny bit longer

plot(c(rep(1,10), rep(2, 10), rep(3,10)),
     c(times1000, times_pareto2, times_pareto3),
     xlab="Pareto tolerance",
     ylab="Run time (mins)",
     main="Pareto tolerance (ranks up to 3)")

plot(c(rep(1,10), rep(2, 10), rep(3,10)),
     c(times1000, times1000neigh2, times1000neigh3),
     xlab="Neighbourhood size",
     ylab="Run time (mins)",
     main="Neighbourhood size (up to 3rd-degree)")

plot(c(rep(1000,10), rep(5000, 10), rep(10000,10)),
     c(times1000, times5000, times10000),
     xlab="Pool size",
     ylab="Run time (mins)",
     main="Pareto tolerance (ranks up to 3)")


final_frontsdf <- rbindlist(final_fronts_educated) %>%
  as.data.frame()
tmp <- final_frontsdf[,grep("site", names(final_frontsdf))]
tmp <- as.data.frame(t(apply(tmp, 1, sort)))
names(tmp) <- paste("site", 1:ncol(tmp), sep="")
final_frontsdf[,grep("site", names(final_frontsdf))] <- tmp
collected <- final_frontsdf %>%
  #dplyr::select(grep("site", names(final_frontsdf), value=TRUE)) %>%
  unique() #%>%

# contour plot but now it's final fronts only .. check this for 10,000 iters
# work out a way to save everything? might just have to be an rds

final_frontsdf <- rbindlist(final_fronts_educated) %>%
  as.data.frame()
agg_pareto <- psel(final_frontsdf, high("sum_pop")*high("sum_risk")) %>%
  arrange(sum_pop)

# a little concerned the same design is in here a bunch of times?
all_sites <- agg_pareto %>%
  dplyr::select(grep("site", names(agg_pareto))) %>%
  unique() %>%
  unlist() %>%
  ftable() %>%
  as.data.frame()
  
agg_map <- vic_objective$potent
values(agg_map)[!is.na(values(agg_map))] <- 0
values(agg_map)[as.numeric(paste(all_sites$.))] <- all_sites$Freq
values(agg_map)[values(agg_map) == 0] <- NA

pot_map <- vic_objective$potent
values(pot_map) <- NA
values(pot_map)[unlist(agg_pareto[1,grep("site", names(agg_pareto))])] <- 1

pot_catch <- vic_objective$potent
values(pot_catch) <- NA
values(pot_catch)[unique(as.vector(catch_membership_mat[unlist(agg_pareto[1,grep("site", names(agg_pareto))]),]))] <- 1

pop_map <- vic_objective$potent
values(pop_map) <- NA
values(pop_map)[unlist(agg_pareto[nrow(agg_pareto), grep("site", names(agg_pareto))])] <- 1

pop_catch <- vic_objective$potent
values(pop_catch) <- NA
values(pop_catch)[unique(as.vector(catch_membership_mat[unlist(agg_pareto[nrow(agg_pareto),grep("site", names(agg_pareto))]),]))] <- 1

mid <- nrow(agg_pareto)/2
mid_map <- vic_objective$potent
values(mid_map) <- NA
values(mid_map)[unlist(agg_pareto[mid,grep("site", names(agg_pareto))])] <- 1

mid_catch <- vic_objective$potent
values(mid_catch) <- NA
values(mid_catch)[unique(as.vector(catch_membership_mat[unlist(agg_pareto[mid,grep("site", names(agg_pareto))]),]))] <- 1


{png("figures/vic_mapped.png",
    width=2500,
    height=2000,
    pointsize=35)
par(mfrow=c(1,2), oma=c(16,0,2,0), mar=c(4,4,2,0))

plot(final_frontsdf$sum_pop, final_frontsdf$sum_risk, 
     xlab="Sum(Human Population Density)",
     ylab="Sum(Potential Risk)",
     cex.lab=1.2, cex.axis=1.2)
lines(rep(max(agg_pareto$sum_pop), 2), c(0, agg_pareto$sum_risk[nrow(agg_pareto)]), col="grey", lty=2, lwd=2)
lines(c(0, agg_pareto$sum_pop[1]), rep(max(agg_pareto$sum_risk), 2), col="grey", lty=2, lwd=2)
points(agg_pareto$sum_pop, agg_pareto$sum_risk, col=brat, pch=16)
lines(agg_pareto$sum_pop, agg_pareto$sum_risk, col=brat, lwd=2)
text(30000, 70, paste("AUF:\n", format(round(pareto_progress_auc(list(agg_pareto)), digits=0), big.mark=",")),
     col=alpha(brat, 0.6), cex=2.5, font=2)
text(existing_hpop, 90, "Existing\n surveillance", col=iddu(2)[2], cex=1.2)
arrows(existing_hpop, 93, existing_hpop, existing_potent-1, col=iddu(2)[2], lwd=2)
points(agg_pareto[c(1,mid,nrow(agg_pareto)), c("sum_pop", "sum_risk")], 
       col=c(berry, "orange", purp), pch=16, cex=1.5)
points(agg_pareto[c(1,mid,nrow(agg_pareto)), c("sum_pop", "sum_risk")], 
       col=c(berry, "orange", purp), cex=3, lwd=2)
points(existing_hpop, existing_potent, col=iddu(2)[2], pch=16, cex=1.5) # make this a little easier to see?
points(existing_hpop, existing_potent, col=iddu(2)[2], cex=3, lwd=2)

par(mar=c(4,2,2,6), bty="n")
plot(agg_map, col=greens(100), axes=FALSE, bty="n", 
     legend.args=list(text="Frequency selected", side=2, line=1, cex=1.2))
plot(st_geometry(vic_shp), add=TRUE)

par(mfrow=c(1,3), oma=c(0,2,50,0), mar=c(0,0,0,0), new=TRUE, mfg=c(1,1))
plot(pot_catch, col=alpha(berry, 0.2), axes=FALSE, bty="n", legend=FALSE) # weird that there are overlapping catchments in here ...
plot(pot_map, add=TRUE, col=berry, legend=FALSE)
# add map of catchment?
plot(st_geometry(vic_shp), add=TRUE)
plot(mid_catch, col=alpha("orange", 0.2), axes=FALSE, bty="n", legend=FALSE)
     #legend.args=list(text="Frequency selected", side=2, line=1))
plot(mid_map, col="orange", add=TRUE, legend=FALSE)
plot(st_geometry(vic_shp), add=TRUE, legend=FALSE)
plot(pop_catch, col=alpha(purp, 0.2), axes=FALSE, bty="n", legend=FALSE)
     #legend.args=list(text="Frequency selected", side=2, line=1))
plot(pop_map, col=purp, add=TRUE, legend=FALSE) # or plot points as in the toy fig ...
plot(st_geometry(vic_shp), add=TRUE)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2,2,2,2), new=TRUE, bty="o", xpd=NA)
empty_plot_for_legend()
subfigure_label(par()$usr, 0,1,"(a)", 1.2)
subfigure_label(par()$usr, 0.55,1,"(b)", 1.2)
subfigure_label(par()$usr, -0.02,0.28,"(c)", 1.2)
subfigure_label(par()$usr, 0.33,0.28,"(d)", 1.2)
subfigure_label(par()$usr, 0.68,0.28,"(e)", 1.2)
dev.off()}














