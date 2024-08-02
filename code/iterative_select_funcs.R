library(rPref)
# write something to extend chain


obj_to_pr = function(obj_out){
  # convert two objective values to probability of accepting
  # assumes we're maximising
  if (obj_out[2] >= obj_out[1]){
    return(1)
  } else {
    return(obj_out[2] / obj_out[1])
  }
}



compare_objectives = function(current, proposed, obj_list, site_catchment_list,
                              dist_mat, raster_stack, pix_weights, find_union){
  # given current and proposed sets of site IDs, report acceptance probability
  #message(current, proposed)
  # first, find union set of pixels in candidate surveillance design
  if (find_union == TRUE){
    current_net_catch = site_catchment_list[unlist(current)] %>%
      unlist() %>%
      unique()
    proposed_net_catch = site_catchment_list[unlist(proposed)] %>%
      unlist() %>%
      unique
  } else {
    current_net_catch = NA
    proposed_net_catch = NA
  }
  
  # for current/proposed, generate list of objective values
  current_objs = lapply(obj_list, function(f){f(current, current_net_catch,
                                                dist_mat, raster_stack,
                                                pix_weights)
  })
  proposed_objs = lapply(obj_list, function(f){f(proposed, proposed_net_catch,
                                                 dist_mat, raster_stack,
                                                 pix_weights)
  })
  # combine objective values
  pr_acc = apply(data.frame(unlist(current_objs),
                         unlist(proposed_objs)),
              1, obj_to_pr) %>%
    prod()
  
  return(list(pr_acc = pr_acc,
              the_details = data.frame(unlist(current_objs),
                                       unlist(proposed_objs))))
}







sim_an_anneal = function(site_catchment_list = list(), 
                          nselect = 0, 
                          niters = 0,
                          raster_stack = stack(),
                          dist_mat = matrix(NA),
                          site_catchment_weights_list = 0,  # hacking this for the poster
                          obj_list = list(),
                          find_union = FALSE,
                          max_temperature = 1,
                         init_sites = c()){
  
  nsites = length(site_catchment_list)
  # initialise object to keep track of selected sites
  outdf = data.frame(matrix(NA, 
                            nrow = niters*nselect + 1, 
                            ncol = nselect))
  if (length(init_sites) == nselect){
    outdf[1, ] = init_sites
  } else {
    outdf[1, ] = sample.int(nsites, 
                            size = nselect,
                            replace = FALSE)
  }
  
  temperatures = seq(1, max_temperature, length.out = nrow(outdf))
  acc_rat = 0
  pr_accs = c()
  deets = data.frame()
  step1 = 0
  step2 = 0
  step3 = 0
  # outer loop: niters (Metropolis)
  for (t in 1: niters){
    # inner loop: nselect (Gibbs)
    for (i in 1: nselect){
      # propose new design
      current = outdf[(t - 1) * nselect + i,]
      proposed = current
      t1 = Sys.time()
      # dropping this step as it's taking up a bunch of time?
      #probs = rep(1, nsites)
      #probs[unlist(current)] = 0
      proposed[i] = sample(x = as.vector(1: nsites)[-unlist(current)],
                           size = 1)#, # maybe prob is taking a while?
                           #prob = sapply(1:nsites, function(x){ifelse(x %in% current, 0, 1)}))
                           #prob = ifelse(1:nsites %in% current, 0, 1),
                           #prob = probs)
      t2 = Sys.time()
      # calculate acceptance probability: obj(current_design, proposed_design)
      pr_acc = compare_objectives(current, proposed, obj_list,
                                  site_catchment_list,
                                  dist_mat, raster_stack, site_catchment_weights_list,
                                  find_union)
      t3 = Sys.time()
      pr_acc$pr_acc = pr_acc$pr_acc ** temperatures[(t - 1) * nselect + i]
      deets = rbind(deets, unlist(pr_acc$the_details))
      pr_acc = pr_acc$pr_acc
      pr_accs = c(pr_accs, unlist(pr_acc))
      t4 = Sys.time()
      
      step1 = step1 + t2 - t1
      step2 = step2 + t3 - t2
      step3 = step3 + t4 - t3
      # accept / reject
      if (runif(1) <= pr_acc){
        outdf[(t - 1) * nselect + i + 1,] = proposed
        acc_rat = acc_rat + 1
      } else {
        outdf[(t - 1) * nselect + i + 1,] = current
      }
    }
  }
  
  # names(deets) = as.vector(outer(c("Current", "Proposed"), 
  #                                1:length(obj_list), paste0))
  names(deets) = c(paste0("Current", 1:length(obj_list)),
                   paste0("Proposed", 1:length(obj_list)))
  
  return(list(outdf=outdf, # contains accepted designs
              acc_rat=acc_rat / (nrow(outdf) - 1),
              pr_accs=pr_accs,
              deets=deets, # contains utilities of current/proposed designs
              step1=step1,step2=step2,step3=step3))
}


# I think this was looking at initial sites, setting a particular proposal 
# distn instead of taking a random sample ...
sim_a_wander = function(site_catchment_list = list(), 
                         nselect = 0, 
                         niters = 0,
                         raster_stack = stack(),
                         connectivity_mat = matrix(NA), # for non-uniform proposals
                         site_catchment_weights_list = list(),
                         obj_list = list(),
                         find_union = FALSE,
                         max_temperature = 1,
                        init_sites = c()){
  
  nsites = length(site_catchment_list)
  # initialise object to keep track of selected sites
  outdf = data.frame(matrix(NA, 
                            nrow = niters*nselect + 1, 
                            ncol = nselect))
  if (length(init_sites) == nselect){
    outdf[1, ] = init_sites
  } else {
    outdf[1, ] = sample.int(nsites, 
                            size = nselect,
                            replace = FALSE)
  }
  temperatures = seq(1, max_temperature, length.out = nrow(outdf))
  acc_rat = 0
  pr_accs = c()
  deets = data.frame()
  # outer loop: niters (Metropolis)
  for (t in 1: niters){
    # inner loop: nselect (Gibbs)
    for (i in 1: nselect){
      # propose new design
      current = outdf[(t - 1) * nselect + i,]
      proposed = current
      #message(setdiff(1: nsites, proposed))
      proposed[i] = sample(x = 1: nsites,
                           size = 1,
                           prob = sapply(1:nsites, function(x){
                             #message(current[i])
                             ifelse(x %in% current, 0, connectivity_mat[x, unlist(current)[i]])}))
      #message(paste(current[i], proposed[i], inv_dist_mat[x,current[i]]))
      
      # calculate acceptance probability: obj(current_design, proposed_design)
      pr_acc = compare_objectives(current, proposed, obj_list,
                                  site_catchment_list,
                                  dist_mat, raster_stack, site_catchment_weights_list,
                                  find_union)
      pr_acc$pr_acc = pr_acc$pr_acc ** temperatures[(t - 1) * nselect + i]
      deets = rbind(deets, unlist(pr_acc$the_details))
      pr_acc = pr_acc$pr_acc
      pr_accs = c(pr_accs, unlist(pr_acc))
      
      # accept / reject
      if (runif(1) <= pr_acc){
        outdf[(t - 1) * nselect + i + 1,] = proposed
        acc_rat = acc_rat + 1
      } else {
        outdf[(t - 1) * nselect + i + 1,] = current
      }
    }
  }
  
  # names(deets) = as.vector(outer(c("Current", "Proposed"), 
  #                                1:length(obj_list), paste0))
  names(deets) = c(paste0("Current", 1:length(obj_list)),
                   paste0("Proposed", 1:length(obj_list)))
  
  return(list(outdf=outdf, # contains accepted designs
              acc_rat=acc_rat / (nrow(outdf) - 1),
              pr_accs=pr_accs,
              deets=deets # contains utilities of current/proposed designs
              
  ))
}


wrap_sim_plots = function(out, 
                          single_site_utility,
                          map_to_single_site_utility = c(),
                          panelled_plot = TRUE, 
                          path = "",
                          main = "",
                          burn_in = 0,
                          trace = FALSE,
                          sandbox = raster(),
                          heat_map_panel = FALSE,
                          sandbox_shp = c(),
                          site_ids = c(),
                          temperature_superimpose = c(), # not yet implemented
                          trace_in_sandbox_panel = FALSE){
  # should work for applications with one objective, 
  # for as many sites per design as I like
  if (path != ""){
    if (trace_in_sandbox_panel == TRUE & length(map_to_single_site_utility) == 0){
      png(path, height=2400, width=2400, pointsize=35)
    } else {png(path, height=1800, width=2400, pointsize=35)}
  }
  
  if (panelled_plot == TRUE){
    if (trace_in_sandbox_panel == TRUE & length(map_to_single_site_utility) == 0){ 
      par(mfrow=c(3, 2), mar=c(4.1,4.1,3.1,4.1)) 
    } else { par(mfrow=c(2, 2), mar=c(4.1,4.1,3.1,4.1)) }
  }

  if (main != ""){par(oma=c(0,0,2.1,0))}
  
  if (burn_in != 0){
    out$outdf = out$outdf[burn_in:nrow(out$outdf),]
    out$deets = out$deets[burn_in:nrow(out$deets),]
  }
  
  nsites = length(single_site_utility)
  
  # partitioned this outside the function ...
  if (length(map_to_single_site_utility) == 0){
    flag = TRUE
    map_to_single_site_utility = sapply(1:nsites, 
                                        function(x){which(order(single_site_utility, 
                                                                decreasing = TRUE) == x)})
    # histogram of sites selected
    out_hist = hist(map_to_single_site_utility[unlist(out$outdf)], breaks=100, 
                    main="Histogram of sites selected", xlab="", xaxt="n")
    
    # superimpose objective values for single sites
    points(map_to_single_site_utility,
           single_site_utility * max(out_hist$counts) / max(single_site_utility), cex=0.8)
    axis(4, at = seq(0, round(max(out_hist$counts), 2), length.out=5), 
         labels = seq(0, round(max(single_site_utility), 2), length.out=5))
    mtext("Site obj", side=4, line = 2, cex=0.8)
  } else { # alternatively, design utility
    flag = FALSE
    hist(out$deets$Current1, breaks=100, 
         main="Histogram of accepted design utilities", xlab="")
  }
  
  
  # acceptance probabilities as simulation goes, with spline over top
  plot(out$pr_accs, cex=0.8,  main="Acceptance probabilities over simulation",
       ylab="Pr(acc)")
  lines(smooth.spline(x=1:length(out$pr_accs), y=out$pr_accs))
  
  # trace plot of current/proposed design utility
  matplot(out$deets[seq(1,nrow(out$deets)),], type="l", 
          main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")
  
  if (trace_in_sandbox_panel == TRUE){
    trace_in_sandbox(out$outdf, sandbox, main="Trace plot in 'geographic' space")
  }
  
  if (heat_map_panel == TRUE){
    # colour pix where there are sites?
    # would be good to colour pix covered by design
    
    # not sure if this crumbles when I don't end up sampling all of the sites
    # works for now ... :/ might need an intermediate step to create empty df then populate
    tmp2 = apply(out$outdf, 2, table)
    tmp2 = rowSums(tmp2)
    sandbox[] = NA
    values(sandbox)[site_ids] = tmp2
    plot(trim(sandbox, padding=5), col=alpha(viridis(100)), main="Frequently selected sites",
         bty="n", xaxt="n", yaxt="n")
    if (length(sandbox_shp) != 0){plot(sandbox_shp, add=TRUE, border=viridis(100)[1])}
  } else {
    # Want to sort below plot by individual site objective
    # matplot(out$outdf, xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?"
    pal = brewer.pal(min(9, ncol(out$outdf)), "Set1")
    if (flag == TRUE){
      plot(0,0, xlab="Iteration", ylab="Site (by single site utility)", 
           pch=16, xlim=c(0,nrow(out$outdf)),
           ylim=c(0, length(single_site_utility)), type="n",
           main = "Which sites get stuck?")
      for (i in 1:ncol(out$outdf)){
        points(1:nrow(out$outdf), map_to_single_site_utility[out$outdf[,i]], 
               col = alpha(pal[i], single_site_utility[out$outdf[,i]] / max(single_site_utility)))
      }
    }
  }
  
  if (main != ""){mtext(main, outer=TRUE)}
  
  if (path != ""){dev.off()}
  
  if (trace == TRUE){plot(out$outdf, type="l")}
}


# atm this grabs coords from rasterToPoints
trace_in_sandbox <- function(outdf, sandbox, 
                             trace_pal=brewer.pal(ncol(outdf), "Set1"),
                             bg_pal=brewer.pal(9,"Purples"),
                             alpha_wt=0.1,
                             path="", main=""){
  if (path != ""){
    png(path, height=1000, width=1000, pointsize=30)
  }
  
  plot(sandbox, col=bg_pal, main=main, legend=FALSE, xaxt="n", yaxt="n", bty="n")
  plot(sandbox, legend.only=TRUE, col=bg_pal, legend.args=list("Pixel obj", side=4, line=2), 
       axis.args=list(at=c(0.5,2.5), labels=c("Low", "High")))
  coords = rasterToPoints(sandbox)
  offsets = seq(-0.01, 0.01, length.out=ncol(outdf))
  lapply(1:ncol(outdf), function(x){
    # kind of love the jittering but it looks ridiculous
    # lines(jitter(coords[outdf[,x], 1:2]), col=trace_pal[x])
    to_plot = jitter(coords[outdf[,x], 1:2])
    # This looks fine ..
    #lines(coords[outdf[,x], 1:2] + offsets[x], col=trace_pal[x])
    
    # Here we go ..    
    sapply(2:nrow(outdf), function(y){
      lines(to_plot[c(y-1,y),] + offsets[x],
            col=alpha(trace_pal[x],
                      alpha=y/nrow(outdf) * (1 - alpha_wt) + alpha_wt), lwd=2)
    })
  })
  
  if (path != ""){dev.off()}
  return(0)
}


post_hoc_pareto <- function(outdf, outobjs, ras, tol=0.2, jitter_fac=1,
                            obj_labs=c()){
  # option for ras to be a stack?
  # offset outobjs by one against outdf
  message(paste0(nrow(unique(outdf)), " unique designs of ", nrow(outdf)))
  # just look at unique designs
  outdf <- unique(outdf)
  outobjs <- outobjs[, seq(1,ncol(outobjs)/2)]
  site_locs <- rasterToPoints(ras)
  
  if (length(obj_labs) == 0){
    obj_labs = paste0("Obj ", 1:ncol(outobjs))
  }
  
  if (is.null(dim(outobjs))){ # one objective
    shortlist_objs <- order(outobjs, decreasing=TRUE)[seq(1,tol * length(outobjs))]
    par(mfrow=c(1,2))
    hist(outobjs, breaks=100, main="Utility of all designs visited", xlab="Utility")
    abline(v=min(outobjs[shortlist_objs]), lty=2, col="red")
    site_polys_in_space_plot(outdf[shortlist_objs,], site_locs, ras, 
                             jitter_fac=jitter_fac, 
                             legend_lower_upper = range(outobjs[shortlist_objs]),
                             panelled = FALSE)
   
  } else { 
    
    p <- high("Current1")
    for (i in 1:(ncol(outobjs) - 1)){
      p <- p * high(outobjs[,i+1])
    }
    #return(p)
    sky1 = psel(outobjs, p, top=nrow(outobjs))
    
    # two objectives
    if (ncol(outobjs) == 2){
      par(mfrow=c(2,2))
      plot(0)
      hist(outobjs[,1], breaks=100)
      hist1 = hist(outobjs[,2], breaks=100, plot=FALSE)
      barplot(tmp$counts, axes=TRUE, space=0, horiz=TRUE)
      plot(sky1[,1:2], xlab=obj_labs[1], ylab=obj_labs[2], 
           col=viridis(max(sky1[,3]))[sky1[,3]]) # add colouring for Pareto
      points(sky1[sky1$.level == 1, 1:2], col="red", pch=16)
      
    }
    
    # more than two objectives
    
    
  }
}

# tmp2=post_hoc_pareto(tmp$outdf, tmp$deets, sandbox, tol=0.01, jitter_fac=0.5)
# columns in outobjs aren't ordered right ..

site_polys_in_space_plot <- function(site_id_df, site_loc_df, ras, jitter_fac,
                                     legend_lower_upper = NA, panelled=TRUE){
  # tolerance applied in outer function
  # probably need to have a play with colours, etc.
  purples = brewer.pal(9, "Purples")
  oranges = brewer.pal(9, "Oranges")
  oranges = colorRampPalette(oranges)(nrow(site_id_df))
  main = ifelse(panelled==TRUE, "", paste0("Best ", nrow(site_id_df), " designs"))
  plot(ras, col=purples, xaxt="n", yaxt="n", legend=FALSE,
       main=main)
  
  nsites <- ncol(site_id_df)
  
  site_id_df <- site_id_df %>%
    rownames_to_column(var="region") %>%
    pivot_longer(!region) %>%
    mutate(lon = jitter(site_loc_df[value, "x"], jitter_fac),
           lat = jitter(site_loc_df[value, "y"], jitter_fac))
    # now for the spatial bit ..
  hulls = site_id_df %>%
    st_as_sf(coords = c("lon", "lat")) %>%
    group_by(region) %>%
    summarize(geometry = st_union(geometry)) %>%
    st_convex_hull()
  
  # plot order reversed so that best designs are on top
  points(site_id_df[nrow(site_id_df):1, c("lon","lat")], col=rep(oranges, each=nsites))
  # something's up with the colouring of the points :/
  plot(hulls[nrow(hulls):1,], col=NA, border=oranges, add=TRUE)
  
  if(length(legend_lower_upper) == 2){ # if upper and lower bounds provided
    make_legend_for_polys(legend_lower_upper, ras, oranges)
  }
  
  # points inside the convex hull don't end up captured by the polygon ... 
  # but I think this is a good thing? So leaving as is.
  
  return(list(df=site_id_df, hulls=hulls))
}

# tmp2 <- site_polys_in_space_plot(tail(tmp$outdf), site_locs, sandbox, jitter_fac=1)

make_legend_for_polys <- function(lower_upper, ras, pal){
  plot(ras, legend.only=TRUE, col=pal, 
       axis.args=list(at=c(minValue(ras), maxValue(ras)),
                      labels=round(lower_upper, 2)),
       legend.args=list("Design utility", side=1, line=2),
       horizontal=TRUE)
  
}

# Genetic_selectic ... attempt 1
# Given dataframe of sites, with lonlats, catchment masks (raster), number of sites,
# number of iterations, one or more objectives (wrt sampling design) - these should be functions!

# obj_to_pr = function(obj_out_curr, obj_out_proposed){
#   if (length(obj_out_curr) != 1 | length(obj_out_proposed) != 1){
#     message("Alert!")
#   }
#   # assumes we're maximising
#   if (obj_out_proposed >= obj_out_curr){
#     return(1)
#   } else {
#     return(obj_out_proposed / obj_out_curr)
#   }
# }
# 
# 
# # don't yet know how to do a vector/list of functions ...
# # each objective takes a DESIGN and returns a POSITIVE NUMBER
# # (unless we want to minimise)
# genetic_selectic1 = function(cands_df, nsites, niters,
#                             obj1 = function(x){return(FALSE)},
#                             obj2 = function(x){return(FALSE)},
#                             obj3 = function(x){return(FALSE)}){
#   # initialise object to keep track of selected sites
#   outdf = data.frame(matrix(NA, 
#                             nrow = niters*nsites + 1, 
#                             ncol = nsites))
#   outdf[1, ] = sample.int(nrow(cands_df), 
#                           size = nsites,
#                           replace = FALSE)
#   acc_rat = 0
#   pr_accs = c()
#   
#   # Put objectives together at this point ...
#   # A bit constrained here on the number of objectives ...
#   nobj = 0
#   if (obj1(1) != FALSE){
#     if (obj2(1) != FALSE){
#       if (obj3(1) != FALSE){
#         compare = function(design1, design2){
#           obj_to_pr(obj1(design1), obj1(design2)) * 
#             obj_to_pr(obj2(design1), obj2(design2)) * 
#             obj_to_pr(obj3(design1), obj3(design2))
#         }
#       } else {
#         obj = function(design1, design2){
#           obj_to_pr(obj1(design1), obj1(design2)) * 
#             obj_to_pr(obj2(design1), obj2(design2))
#         }
#       } 
#     } else {
#       obj = function(design1, design2){
#         obj_to_pr(obj1(design1), obj1(design2))
#       }
#     } 
#   } else {
#     message("False start")
#     return(FALSE)
#   }
#   
#   # outer loop: nits (Metropolis)
#   for (t in 1: niters){
#     # inner loop: nsites (Gibbs)
#     for (i in 1: nsites){
#       # propose new design
#       tmp = outdf[(t - 1)*nsites+ i,]
#       tmp[i] = sample(setdiff(1:nrow(cands_df), tmp), # could use probs argument instead of setdiff
#                       size = 1)
#       # calculate acceptance probability: obj(current_design, proposed_design)
#       return(cands_df[unlist(outdf[(t - 1)*nsites + i,]),])
#       pr_acc = obj(cands_df[outdf[(t - 1)*nsites + i,],], 
#                    cands_df[tmp,])
#       pr_accs = c(pr_accs, unlist(pr_acc))
#       # accept / reject
#       if (runif(1) <= pr_acc){
#         outdf[(t - 1)*nsites + i + 1,] = tmp
#         acc_rat = acc_rat + 1
#       } else {
#         outdf[(t - 1)*nsites + i + 1,] = outdf[(t - 1)*nsites + i,]
#       }
#     }
#   }
#   return(list(outdf=outdf,
#               acc_rat=acc_rat / (nrow(outdf) - 1),
#               pr_accs=pr_accs))
# }
# 
# indf = data.frame(matrix(rnorm(15), nrow=5, ncol=3))
# indf[4,] = rnorm(3, mean=3)
# 
# tmp = genetic_selectic(indf,
#                  3, 5, 
#                  obj1 = sum, 
#                  obj2 = function(x){if (length(x) > 1){x[3]} else {TRUE}})