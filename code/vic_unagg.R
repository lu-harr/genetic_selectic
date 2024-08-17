vic_shp <- states %>%
  filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)

# need 20km to do third-degree neighbours ...
vic_shadow <- vic_shp %>%
  st_buffer(dist = 25000) %>% # 80% certain we're in metres
  st_simplify(dTolerance = 1000)

vic_objective <- stack(raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif'),
                       raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif')) %>%
  # aggregate(AGG_FACTOR) %>%
  crop(extent(vic_shadow)) %>%
  mask(vic_shadow)
names(vic_objective) <- c("potent", "hpop")

tmp <- vic_objective
vic_objective <- mask(vic_objective, vic_shp)
vic_objective$buffer_potent <- tmp$potent
vic_objective$buffer_hpop <- tmp$hpop

id_ras <- vic_objective$potent
id_ras[] <- 1:ncell(id_ras)

vic_mozzies <- subset(all_mozzies, data_source == "Vic") %>%
  filter(as.Date(date) > as.Date("2022-07-01")) %>% 
  # consider only 2022-2023 trapping season
  group_by(longitude, latitude) %>%
  summarise(ntimes = n()) %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise()

neigh_mat <- focalWeight(id_ras, 0.05, "circle") # 5kms ish ...
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

# 789,756 total sites
# 327,652 vic sites
# 244 rows in vic_mozzies

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

nselect <- 200

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

plot(vic_objective$potent)
for (i in 1:nrow(educated_guess)){
  tmp <- which(vic_sites$id %in% educated_guess[i,])
  points(vic_sites[tmp, c("x","y")])
}
# this now takes considerably longer (probably the which())

niters = 100
{set.seed(834903)
  tstart1 <- Sys.time()
  tmp <- genetic_algot(site_ids = vic_sites$id,
                       nselect = nselect, 
                       poolsize = 1000,
                       niters = 100,
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
  tend1 <- Sys.time()} #  mins to do 100 iters
tend1 - tstart1
# 29 minutes
# so this takes 17 times longer to do ..
# so 10,000 iters would take 42.5 hours
# which is worth doing but not feasible to do on my lappy








