library(raster)
library(RColorBrewer)
library(scales)
library(dplyr)
library(viridisLite)

plotpath = "~/Desktop/knowlesi/multi_site/exploratory_plots/"


##############################################################################
# TAKE TWO ...
# lets set up a toy problem that makes sense and allows me to try components individually
sandbox = raster(matrix(runif(256), nrow=16, ncol=16))
sandbox[3:6, 3:6] = 2
sandbox[7:10,7:10] = 3
tmp = sandbox
values(tmp) = 1:ncell(sandbox)
sandcastle = stack(sandbox,
                   exp(sandbox),
                   tmp)
names(sandcastle) = c("one", "two", "id")

# playsites = c(25, 69, 110, 200, 135, 74, 208)
# # these represent pixel ids .. obvs can automate this :)
# playcatch = list(c(9,10,11,24,25,26,40,41,42),
#                  c(53,54,55,68,69,70),
#                  c(93,94,95,109,110,111,125,126,127),
#                  c(184,183,185,199,200,201,215,216,217),
#                  c(135,136,137, 119, 121, 120, 151, 153, 152),
#                  c(73,74,75,58,57,59,42,41,43),
#                  c(208,207,206,192,191,190,176,175,174))
plot(sandbox)
# leave distmat for now
# there's got to be a way to do this using focal() or something similar

playsites = 1:ncell(sandbox)
# doesn't include itself :/
playcatch = lapply(1:ncell(sandbox), function(x){adjacent(sandbox, x, directions=8, 
                                                          pairs=FALSE, include=TRUE)})



# for example ...
# obj_list[[1]](pix_in_catch = c(53,54,55,68,69,70), raster_stack = sandcastle)
# 
# sapply(1:length(playcatch), function(x){obj_list[[1]](pix_in_catch = playcatch[[x]],
#                                                       raster_stack = sandcastle)})
# 
# # also works over whole list ..
# lapply(obj_list, function(f){obj_to_pr(f(pix_in_catch = c(53,54,55,68,69,70),
#                                          raster_stack = sandcastle))}
# )

# Here's an example objective_list
objective_list = list(
  obj1 = function(site_ids = c(), 
                  pix_in_catch = c(), # a vector of all pixels across all catchments in design
                  dist_mat = matrix(NA), 
                  raster_stack = stack(), 
                  pix_weights = c()){
    # sum objective values of pix in catch (first stack in raster)
    sum(raster_stack[[1]][pix_in_catch])
   }#,
  #  obj2 = function(site_ids = c(), pix_in_catch = c(),
  #                 dist_mat = matrix(NA),
  #                 raster_stack = stack(),
  #                 pix_weights = c()){
  #   message(pix_in_catch)
  #   mean(raster_stack[[1]][pix_in_catch])
  # },
  # obj3 = function(site_ids = c(), pix_in_catch = c(),
  #                 dist_mat = matrix(NA),
  #                 raster_stack = stack(),
  #                 pix_weights = c()){
  #   values(raster_stack[[1]])[site_ids]
  # }
)


plot(sandbox, col=brewer.pal(9, "Purples"))
single_site_utility = sapply(1:length(playcatch), 
                             function(x){objective_list[[1]](pix_in_catch = playcatch[[x]],
                                                             raster_stack = sandcastle)})
points(rasterToPoints(sandbox)[,1:2],
       col=alpha("orange", single_site_utility/max(single_site_utility)),
       pch=16)



tmp = sim_an_anneal(site_catchment_list = playcatch,
                    nselect = 3, 
                    niters = 500, 
                    raster_stack = sandcastle,
                    obj_list = objective_list,
                    find_union = TRUE,
                    max_temperature = 1)

wrap_sim_plots(tmp, two_by_two = TRUE, single_site_utility = single_site_utility,
               main="Select 3 sites from sandbox: no annealing",
               path = paste0(plotpath,"sim_max_temp1.png"))

# add to wrap sim - site *combinations*

# increasing temperature does slow things down ..
tmp = sim_an_anneal(site_catchment_list = playcatch,
                    nselect = 3, 
                    niters = 500, 
                    raster_stack = sandcastle,
                    obj_list = objective_list,
                    find_union = TRUE,
                    max_temperature = 10)

wrap_sim_plots(tmp, two_by_two = TRUE,  single_site_utility = single_site_utility,
               main="Select 3 sites from sandbox: with annealing (temp 1 -> 10)",
               path = paste0(plotpath,"sim_max_temp10.png"))


# A couple of examples with max temp ramped up to 50 ..
set.seed(834903)
tmp = sim_an_anneal(site_catchment_list = playcatch,
                    nselect = 3, 
                    niters = 500, 
                    raster_stack = sandcastle,
                    obj_list = objective_list,
                    find_union = TRUE,
                    max_temperature = 50)
set.seed(3)
tmp2 = sim_an_anneal(site_catchment_list = playcatch,
                    nselect = 3, 
                    niters = 500, 
                    raster_stack = sandcastle,
                    obj_list = objective_list,
                    find_union = TRUE,
                    max_temperature = 50)
# set.seed(7)
# tmp3 = sim_an_anneal(site_catchment_list = playcatch,
#                      nselect = 3, 
#                      niters = 500, 
#                      raster_stack = sandcastle,
#                      obj_list = objective_list,
#                      find_union = TRUE,
#                      max_temperature = 50)

png(paste0(plotpath, "sim_maxtemp50_compare.png"),
    height=1800, width=2400, pointsize=35)
par(mfrow=c(2,2))

matplot(tmp$deets[seq(1,nrow(tmp$deets)),], type="l",
        main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")

pal = brewer.pal(ncol(tmp$outdf), "Set1")
plot(0,0, xlab="Iteration", ylab="Site (by single site utility)", 
     pch=16, xlim=c(0,nrow(tmp$outdf)),
     ylim=c(0, length(site_objs)), type="n",
     main = "Which sites get stuck?")
for (i in 1:ncol(tmp$outdf)){
  points(1:nrow(tmp$outdf), map_to_utility[tmp$outdf[,i]], 
         col = alpha(pal[i], site_objs[tmp$outdf[,i]] / max(site_objs)))
}
#matplot(tmp$outdf, xlab="Iteration", ylab="Site ID", pch=1, 
#        main="Which sites get stuck?")

matplot(tmp2$deets[seq(1,nrow(tmp2$deets)),], type="l",
        main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")

plot(0,0, xlab="Iteration", ylab="Site (by single site utility)", 
     pch=16, xlim=c(0,nrow(tmp2$outdf)),
     ylim=c(0, length(site_objs)), type="n",
     main = "Which sites get stuck?")
for (i in 1:ncol(tmp2$outdf)){
  points(1:nrow(tmp2$outdf), map_to_utility[tmp2$outdf[,i]], 
         col = alpha(pal[i], site_objs[tmp2$outdf[,i]] / max(site_objs)))
}

# matplot(tmp3$deets[seq(1,nrow(tmp3$deets)),], type="l",
#         main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")
# 
# plot(0,0, xlab="Iteration", ylab="Site (by single site utility)", 
#      pch=16, xlim=c(0,nrow(tmp3$outdf)),
#      ylim=c(0, length(site_objs)), type="n",
#      main = "Which sites get stuck?")
# for (i in 1:ncol(tmp3$outdf)){
#   points(1:nrow(tmp3$outdf), map_to_utility[tmp3$outdf[,i]], 
#          col = alpha(pal[i], site_objs[tmp3$outdf[,i]] / max(site_objs)))
# }
#matplot(tmp2$outdf, xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?")
dev.off()





















# The below code is now wrapped up in wrap_sim_plots():
# {par(mfrow=c(2, 2))
# 
# out_hist = hist(unlist(tmp$outdf), breaks=100, main="Histogram of sites selected", xlab="Site ID")
# site_objs = sapply(1:length(playcatch), 
#                    function(x){objective_list[[1]](pix_in_catch = playcatch[[x]],
#                                                    raster_stack = sandcastle)})
# points(site_objs * max(out_hist$counts) / max(site_objs), cex=0.8)
# axis(4, at = seq(0, max(out_hist$counts), length.out=5), 
#      labels = seq(0, max(site_objs), length.out=5))
# mtext("Site obj", side=4, line = 2, cex=0.8)
# 
# 
# #hist(tmp$pr_accs, breaks=50, main="Histogram of acceptance probabilities", xlab="Pr(accept)")
# plot(tmp$pr_accs, cex=0.8,  main="Acceptance probabilities over simulation",
#      ylab="Pr(acc)")
# lines(smooth.spline(x=1:length(tmp$pr_accs), y=tmp$pr_accs))
# 
# matplot(tmp$deets[seq(1,nrow(tmp$deets)),]/37, type="l", 
#         main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")
# 
# # Want to sort below plot by individual site objective
# matplot(tmp$outdf, xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?")}
# 
# plot(seq(1, 10, length.out = nrow(tmp$outdf)-1), tmp$pr_accs)



# Do some reading about acceptance probabilities (Muller papers)
# Do that before you do anything else!!!!

################################################################################
# Enumerate all combinations and check overlaps are not advantageous..
# surface of x/y sorted by single site utility, coloured by joint utility

# site_catchment_list = playcatch,
# nselect = 3, 
# raster_stack = sandcastle,
# obj_list = objective_list,
# find_union = TRUE,
# max_temperature = 50

#playcatch # list of catchment ownership

#single_site_utility = sapply(1:length(playcatch), function(x){objective_list[[1]](pix_in_catch = playcatch[[x]],
#                                                                                  raster_stack = sandcastle)})
#map_to_utility = sapply(1:length(playsites), function(x){which(order(single_site_utility, decreasing = TRUE) == x)})
# matplot(tmp$outdf[order(single_site_utility, decreasing = TRUE)[1:20],], 
#         xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?")

pairwise_objs = sapply(1:length(playcatch), function(row){
  sapply(1:length(playcatch), function(col){
    return(objective_list[[1]](site_ids = c(row, col),
                  pix_in_catch = unique(c(playcatch[[row]], playcatch[[col]])),
                  raster_stack = sandcastle))
  })
})

diag(pairwise_objs) = NA

# now order rows and columns by site utility?
ordered_pairwise_objs = pairwise_objs[order(single_site_utility, decreasing = TRUE),]
ordered_pairwise_objs = ordered_pairwise_objs[,order(single_site_utility, decreasing = TRUE)]

neigh_mat = sapply(1:length(playcatch), function(row){
  return(ifelse(1:length(playcatch) %in% playcatch[[row]], 1, 0))
})
diag(neigh_mat) = NA

zoomed_pairs = ordered_pairwise_objs[1:64,1:64]
plot(raster(zoomed_pairs))

{png(paste0(plotpath, "sandbox_2combs_surface.png"),
    height=1500, width=1500, pointsize=35)
par(mfrow=c(2,2), mar=c(2,2,2,2), bty="n")
plot(raster(pairwise_objs), col=viridis(100), xaxt="n", yaxt="n", 
     main="Joint utility (combinations of 2) by site ID")

plot(raster(ordered_pairwise_objs), col=viridis(100), xaxt="n", yaxt="n",
     main="Joint utility (combs of 2) by indiv site utility")
lines(c(0,0.25,0.25,0,0), c(0.75,0.75,1,1,0.75), col="red", lty=2, lwd=4)

plot(raster(zoomed_pairs),
     col=viridis(100), xaxt="n", yaxt="n",
     main="Zoom on best sites in middle panel")

# panel highlighting combs of adjacent sites?
plot(raster(neigh_mat[1:64, 1:64]), xaxt="n", yaxt="n", main="Neighbours")

dev.off()}

# margin plot of site obj so we can see the 24s next to everything else?

# Need to remove diagonal, add site back into playcatch

# How about larger catchments? Does overlap start to become important?



