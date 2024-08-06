vic_shp <- states %>%
  filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)

vic_shadow <- vic_shp %>%
  st_buffer(dist = 5000) %>% # 80% certain we're in metres
  st_simplify(dTolerance = 5000)
# ugly ! ah well !

vic_objective <- stack(raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif'),
                       raster('~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif')) %>%
                         aggregate(AGG_FACTOR) %>%
                         crop(extent(vic_shp)) %>%
                         mask(vic_shadow)
names(vic_objective) <- c("potent", "hpop")
tmp <- vic_objective$potent
vic_objective <- mask(vic_objective, vic_shp)
vic_objective$potent_buffer <- tmp

id_ras <- vic_objective$potent
id_ras[] <- 1:ncell(id_ras)

# I guess compare performances by plotting fronts next to each other?
# Stop when there hasn't been a change for 10 iters?

# try unaggregated case first ...

# raster level assessment of existing trapping!
# all_mozzies is read in in sa4_appraisal ... should transport into main.R

# just want unique locations
vic_mozzies <- subset(all_mozzies, data_source == "Vic") %>%
  group_by(longitude, latitude) %>%
  summarise(ntimes = n()) %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise()
# now need to re-assign latlons ... 

  
raster::extract(id_ras, vic_mozzies[,c("longitude", "latitude")])
  


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