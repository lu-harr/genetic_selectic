# this feels silly but I'm just going to repeat everything I have in the vic script for wa
# try doing some preliminary runs on spartan ... I don't have the heart to run this all on my lappy

wa_shp <- states %>%
  filter(STE_NAME21 == "Western Australia") %>%
  st_simplify(dTolerance = 1000)

# need 20km to do third-degree neighbours ...
wa_shadow <- wa_shp %>%
  st_buffer(dist = 25000) %>% # 80% certain we're in metres
  st_simplify(dTolerance = 1000)
# ugly ! ah well !

wa_objective <- objective_rasters %>%
  aggregate(AGG_FACTOR) %>%
  crop(extent(wa_shadow)) %>%
  mask(wa_shadow)
names(wa_objective) <- c("potent", "hpop")

tmp <- wa_objective
wa_objective <- mask(wa_objective, wa_shp)
wa_objective$buffer_potent <- tmp$potent
wa_objective$buffer_hpop <- tmp$hpop

id_ras <- wa_objective$potent
id_ras[] <- 1:ncell(id_ras)

# check all_mozzies is what I think it is ...
# wa_mozzies <- subset(all_mozzies, data_source == "WA") %>%
#   filter(as.Date(date) > as.Date("2022-07-01")) %>% 
#   # consider only 2022-2023 trapping season
#   group_by(longitude, latitude) %>%
#   summarise(ntimes = n()) %>%
#   mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
#   group_by(pix) %>%
#   summarise()

neigh_mat <- focalWeight(id_ras, 0.12, "circle") # queens case
neigh_mat[!neigh_mat == 0] = 1
catchment_stack <- terra::focal(terra::rast(id_ras), neigh_mat, fun=c)
catchment_stack <- subset(catchment_stack, which(neigh_mat != 0))
catch_membership_mat <- values(catchment_stack, mat=TRUE)

site_ids <- as.data.frame(rasterToPoints(id_ras))
site_ids$potent <- values(wa_objective$buffer_potent)
site_ids$hpop <- values(wa_objective$buffer_hpop)
site_ids$id <- 1:ncell(wa_objective)

wa_sites <- site_ids[!is.na(values(wa_objective$potent)),]
nrow(wa_sites)


# objective_func = function(x, catch_mem, vec){
#   sum(vec[unique(as.vector(catch_mem[unlist(x),]))], na.rm=TRUE)
# }
# 
# existing_potent <- objective_func(wa_mozzies$pix, 
#                                   catch_membership_mat, 
#                                   wa_objective$buffer_potent)
# existing_hpop <- objective_func(wa_mozzies$pix,
#                                 catch_membership_mat,
#                                 wa_objective$buffer_hpop)


################################################################################
# let's deploy GA 

nselect <- 100 # leave the same as Vic even though it's a much larger area ........

set.seed(834904)
starting_point <- matrix(sample(wa_sites$id, 100000, replace=TRUE), ncol=nselect) %>%
  as.data.frame()
names(starting_point) <- paste("site", 1:nselect, sep="")

# educated guess starting point
educated_guess <- matrix(c(sample(wa_sites$id[order(wa_sites$potent, 
                                                     decreasing = TRUE)][1:200], 50000, replace = TRUE),
                           sample(wa_sites$id[order(wa_sites$hpop, 
                                                     decreasing = TRUE)][1:200], 50000, replace = TRUE)),
                         ncol=nselect, byrow=TRUE) %>%
  as.data.frame()
names(educated_guess) <- paste("site", 1:nselect, sep="")

# greedy_sites <- educated_guess %>%
#   unlist() %>%
#   ftable() %>%
#   as.data.frame()
# 
# wa_greedy_map <- wa_objective$potent
# values(wa_greedy_map)[!is.na(values(wa_greedy_map))] <- 0
# values(wa_greedy_map)[as.numeric(paste(greedy_sites$.))] <- greedy_sites$Freq
# wa_greedy_map <- trim(wa_greedy_map)
# values(wa_greedy_map)[values(wa_greedy_map) == 0] <- NA
# 
# # I guess this means I now need to whack this in with the victorian results
# message("here's a greedy start figure")
# {png("figures/wa_greedy_start.png",
#      height=1200,
#      width=1700,
#      pointsize=40)
#   par(mar=c(2,2,4.1,4.1), bty="n", xpd=NA)
#   plot(wa_greedy_map, col=greens(100), axes=FALSE, bty="n",
#        legend.args=list(text="Frequency", 2, line=1),
#        main="Greedy starting pool")
#   par(xpd=NA)
#   plot(st_geometry(wa_shp), add=TRUE)
#   dev.off()}


################################################################################
# baseline run: poolsize 1000







