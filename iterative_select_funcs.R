

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
                           prob = sapply(1:nsites, function(x){ifelse(x %in% current, 0, 1)}))
      
      # calculate acceptance probability: obj(current_design, proposed_design)
      pr_acc = compare_objectives(current, proposed, obj_list,
                                  site_catchment_list,
                                  dist_mat, raster_stack, site_catchment_weights_list,
                                  find_union)
      pr_acc$pr_acc = pr_acc$pr_acc ** temperatures[(t - 1) * nselect + i]
      deets = rbind(deets, pr_acc$the_details)
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
  
  return(list(outdf=outdf, # contains accepted designs
              acc_rat=acc_rat / (nrow(outdf) - 1),
              pr_accs=pr_accs,
              deets=deets # contains utilities of current/proposed designs
              ))
}


sim_a_wander = function(site_catchment_list = list(), 
                         nselect = 0, 
                         niters = 0,
                         raster_stack = stack(),
                         dist_mat = matrix(NA),
                         site_catchment_weights_list = list(),
                         obj_list = list(),
                         find_union = FALSE,
                         max_temperature = 1,
                        init_sites = c(),
                        proposal_sd = 5){
  
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
                             ifelse(x %in% current, 0, dnorm(x, as.numeric(current[i]), sd = proposal_sd))}))
      
      # calculate acceptance probability: obj(current_design, proposed_design)
      pr_acc = compare_objectives(current, proposed, obj_list,
                                  site_catchment_list,
                                  dist_mat, raster_stack, site_catchment_weights_list,
                                  find_union)
      pr_acc$pr_acc = pr_acc$pr_acc ** temperatures[(t - 1) * nselect + i]
      deets = rbind(deets, pr_acc$the_details)
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
  
  return(list(outdf=outdf, # contains accepted designs
              acc_rat=acc_rat / (nrow(outdf) - 1),
              pr_accs=pr_accs,
              deets=deets # contains utilities of current/proposed designs
  ))
}




wrap_sim_plots = function(out, 
                          single_site_utility,
                          nsites = 256,
                          two_by_two = TRUE, 
                          path = "",
                          trace = FALSE,
                          temperature_superimpose = c(),
                          main = ""){
  # should work for applications with one objective, 
  # for as many sites per design as I like
  if (path != ""){png(path, height=1800, width=2400, pointsize=35)}
  
  if (two_by_two == TRUE){par(mfrow=c(2, 2), mar=c(4.1,4.1,3.1,4.1))}
  
  if (main != ""){par(oma=c(0,0,2.1,0))}
  
  #nsites = max(out$outdf) #for some reason this isn't right
  map_to_single_site_utility = sapply(1:nsites, 
                                      function(x){which(order(single_site_utility, 
                                                              decreasing = TRUE) == x)})
  
  # histogram of sites selected
  out_hist = hist(map_to_single_site_utility[unlist(out$outdf)], breaks=100, 
                  main="Histogram of sites selected", xlab="", xaxt="n")
  
  # superimpose objective values for single sites
  points(map_to_single_site_utility,
         single_site_utility * max(out_hist$counts) / max(site_objs), cex=0.8)
  axis(4, at = seq(0, max(out_hist$counts), length.out=5), 
       labels = seq(0, max(single_site_utility), length.out=5))
  mtext("Site obj", side=4, line = 2, cex=0.8)
  
  # acceptance probabilities as simulation goes, with spline over top
  plot(out$pr_accs, cex=0.8,  main="Acceptance probabilities over simulation",
       ylab="Pr(acc)")
  lines(smooth.spline(x=1:length(out$pr_accs), y=out$pr_accs))
  
  # trace plot of current/proposed design utility
  matplot(out$deets[seq(1,nrow(out$deets)),], type="l", 
          main="Toy problem: current/proposed designs", xlab="Iteration", ylab="Utility")
  
  # Want to sort below plot by individual site objective
  # matplot(out$outdf, xlab="Iteration", ylab="Site ID", pch=1, main="Which sites get stuck?"
  pal = brewer.pal(ncol(out$outdf), "Set1")
  plot(0,0, xlab="Iteration", ylab="Site (by single site utility)", 
       pch=16, xlim=c(0,nrow(out$outdf)),
       ylim=c(0, length(site_objs)), type="n",
      main = "Which sites get stuck?")
  for (i in 1:ncol(out$outdf)){
    points(1:nrow(out$outdf), map_to_single_site_utility[out$outdf[,i]], 
           col = alpha(pal[i], site_objs[out$outdf[,i]] / max(site_objs)))
  }
  
  if (main != ""){mtext(main, outer=TRUE)}
  
  if (path != ""){dev.off()}
  
  if (trace == TRUE){plot(out$outdf, type="l")}
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