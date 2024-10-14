# set myself up for Victorian runs
# (watch out there was some masking happening here from other loaded libraries :( )
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
nruns <- 10
starting_pool_size <- 1000

set.seed(834904)
# starting_point <- matrix(sample(vic_sites$id, 100000, replace=TRUE), ncol=nselect) %>%
#   as.data.frame()
# names(starting_point) <- paste("site", 1:nselect, sep="")

starting_point <- lapply(1:nruns, function(x){
  matrix(sample(vic_sites$id, nselect*starting_pool_size, replace=TRUE),
         ncol = nselect) %>%
    as.data.frame() %>%
    setNames(paste("site", 1:nselect, sep=""))
})

educated_guess <- matrix(c(sample(vic_sites$id[order(vic_sites$potent, 
                                                     decreasing = TRUE)][1:200], 
                                  starting_pool_size*nselect/2, replace = TRUE),
                           sample(vic_sites$id[order(vic_sites$hpop, 
                                                     decreasing = TRUE)][1:200], 
                                  starting_pool_size*nselect/2, replace = TRUE)),
                         ncol=nselect, byrow=TRUE) %>%
  as.data.frame()
names(educated_guess) <- paste("site", 1:nselect, sep="")

naive_pop <- vic_sites$id[order(vic_sites$hpop, decreasing = TRUE)][1:100]
naive_risk <- vic_sites$id[order(vic_sites$potent, decreasing = TRUE)][1:100]

# greedy_sites <- educated_guess %>%
#   unlist() %>%
#   ftable() %>%
#   as.data.frame()
# 
# greedy_map <- vic_objective$potent
# values(greedy_map)[!is.na(values(greedy_map))] <- 0
# values(greedy_map)[as.numeric(paste(greedy_sites$.))] <- greedy_sites$Freq
# greedy_map <- trim(greedy_map)
# values(greedy_map)[values(greedy_map) == 0] <- NA
# 
# leg_minmax <- range(values(greedy_map), values(wa_greedy_map), na.rm=TRUE)
# 
# {png("figures/greedy_start.png",
#      height=1500,
#      width=2200,
#      pointsize=40)
#   par(mfrow=c(1,2), mar=c(0,0.2,3.1,1.1), xpd=NA,
#       oma=c(0,2,0,4), bty="n")
#   plot(greedy_map, col=greens(100), axes=FALSE,
#        breaks=seq(leg_minmax[1], leg_minmax[2], length.out=101),
#        legend=FALSE, legend.mar=0)
#   par(xpd=NA)
#   plot(st_geometry(vic_shp), add=TRUE)
#   #par(oma=c(0,0,0,0), new=TRUE)
#   plot(wa_greedy_map, col=greens(100), axes=FALSE,
#        breaks=seq(leg_minmax[1], leg_minmax[2], length.out=101),
#        legend=FALSE)
#   plot(st_geometry(wa_shp), add=TRUE)
#   par(mfrow=c(1,1), new=TRUE, mar=c(0,4.1,3.1,2.1), oma=c(0,0,0,2))
#   plot(wa_greedy_map, col=greens(100), legend.only=TRUE,
#        legend.args=list(text="Frequency", 2, line=1),
#        breaks=seq(leg_minmax[1], leg_minmax[2], length.out=101),
#        axis.args=list(at=seq(250,500,50), labels=seq(250,500,50)),
#        legend.mar=2)
#   mtext("Greedy starting pools", 3, font=2, cex=1.4, line=1)
#   par(new=TRUE, xpd=NA, bty="n")
#   empty_plot_for_legend()
#   subfigure_label(par()$usr, -0.05, 0.95, "(a)")
#   subfigure_label(par()$usr, 0.5, 0.95, "(b)")
#   dev.off()}



