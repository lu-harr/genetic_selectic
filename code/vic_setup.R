# set myself up for Victorian runs

vic_shp <- states %>%
  filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)

# need 20km to do third-degree neighbours ...
vic_shadow <- vic_shp %>%
  st_buffer(dist = 25000) %>% # 80% certain we're in metres
  st_simplify(dTolerance = 1000)
# ugly ! ah well !

vic_objective <- objective_rasters %>%
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