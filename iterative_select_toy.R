

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
playcatch = lapply(1:ncell(sandbox), function(x){adjacent(sandbox, x, directions=8, pairs=FALSE)})

sandbox[unlist(playcatch)] = 4
plot(sandbox)

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
  obj1 = function(site_ids = c(), pix_in_catch = c(), 
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
  # }
)

tmp = sim_an_anneal(site_catchment_list = playcatch,
                    nselect = 2, 
                    niters = 500, 
                    raster_stack = sandcastle,
                    obj_list = objective_list,
                    find_union = TRUE,
                    max_temperature = 20)



{par(mfrow=c(2, 2))

out_hist = hist(unlist(tmp$outdf), breaks=100, main="Histogram of sites selected", xlab="Site ID")
site_objs = sapply(1:length(playcatch), 
                   function(x){objective_list[[1]](pix_in_catch = playcatch[[x]],
                                                   raster_stack = sandcastle)})
points(site_objs * max(out_hist$counts) / max(site_objs), cex=0.8)
axis(4, at = seq(0, max(out_hist$counts), length.out=5), 
     labels = seq(0, max(site_objs), length.out=5))
mtext("Site obj", side=4, line = 2, cex=0.8)


#hist(tmp$pr_accs, breaks=50, main="Histogram of acceptance probabilities", xlab="Pr(accept)")
plot(tmp$pr_accs, cex=0.8,  main="Acceptance probabilities over simulation",
     ylab="Pr(acc)")
lines(smooth.spline(x=1:length(tmp$pr_accs), y=tmp$pr_accs))

matplot(tmp$deets[seq(1,nrow(tmp$deets)),]/37, type="l", 
        main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")

# Want to sort below plot by individual site objective
matplot(tmp$outdf, xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?")}

plot(seq(1, 10, length.out = nrow(tmp$outdf)-1), tmp$pr_accs)



# Do some reading about acceptance probabilities (Muller papers)
# Do that before you do anything else!!!!
single_site_utility = sapply(1:length(playcatch), function(x){objective_list[[1]](pix_in_catch = playcatch[[x]],
                                                                                  raster_stack = sandcastle)})
map_to_utility = sapply(1:length(playsites), function(x){which(order(single_site_utility) == x)})
matplot(tmp$outdf[order(single_site_utility),], xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?")






