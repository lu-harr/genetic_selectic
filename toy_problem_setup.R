# toy problem set up ....
AGG_FACTOR = 10

guelphia_potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif') %>%
  crop(c(141, 146, -37.9, -32.9)) %>%
  aggregate(AGG_FACTOR) %>%
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
enumerated <- combn(ncell(toy_objective$potent), 5)

tmp <- apply(enumerated, 2, function(x){
    # can't count sites in more than one catchment twice:
    sum(site_ids$hpop[unique(as.vector(catch_membership_mat[x,]))], na.rm=TRUE)
  })

write.csv(cbind(enumerated, tmp),
          "~/Desktop/knowlesi/multi_site/output/toy_enumerated.csv")

tmp2 <- apply(enumerated, 2, function(x){
  sum(site_ids$potent[unique(as.vector(catch_membership_mat[x,]))], na.rm=TRUE)
})

enumerated <- combn(10, 5, function(x){
  c(x, sum(x), sd(x))
})

# illustrate that genetic algot finds Pareto front
