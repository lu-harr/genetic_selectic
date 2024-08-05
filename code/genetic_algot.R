# genetic algot functions in here ... do setup + fenangling elsewhere pls

genetic_algot <- function(site_ids, objective_list, nselect, 
                           potential_vec,
                           pop_vec,
                           sandpit,
                           catchment_matrix=matrix(0),
                           neighbourhood_matrix=matrix(0),
                           pool=matrix(0), 
                           niters=10, 
                           poolsize=20, 
                           sample_method="random",
                           box_extent = c(2, 28, 2, 9.66),
                           top_level = 1){
  # site_ids is the vector of all site IDs (from the raster, indices of catch_potential_vec and distmat)
  # catch_potential_vec is vector of *catchment* potential risk
  # neighbourhood_matrix is used to control addition of new designs to pool via "neighbours" method
  # pool is matrix or dataframe of sites only - we add objective values in a sec
  sandpit[!is.na(sandpit)] <- 0
  sandcastle <- stack()
  
  # calculate objective for things in initial pool?
  if (ncol(pool) == nselect){
    obj <- data.frame(sum_risk = apply(pool, 1, function(x){objective_list[[1]](site_ids=unique(as.vector(catchment_matrix[unlist(x),])), 
                                                                                vector_to_be_indexed=potential_vec)}),
                      sum_pop = apply(pool, 1, function(x){objective_list[[1]](site_ids=unique(as.vector(catchment_matrix[unlist(x),])),
                                                                               vector_to_be_indexed=pop_vec)}),
                      addition = rep(FALSE, nrow(pool)))
    pool <- cbind(obj, pool)
  }
  
  # keep track of addition and removal from pool
  progress <- data.frame(new_additions_to_pool=rep(NA,niters),
                         total_on_pareto_front=rep(NA,niters),
                         rejects=rep(NA,niters))
  png("~/Desktop/knowlesi/multi_site/output/toy_movie/iter0.png", width=2200, height=1400, pointsize=40)
  par(mfrow=c(1,2), mar=c(5.1,4.1,2,2), oma=c(0,0,4.1,1.1))
  plot(pool$sum_risk, log(pool$sum_pop), pch=16, col=apple,
       xlim=box_extent[1:2], ylim=box_extent[3:4],
       xlab="Total potential risk captured", ylab="Total human population captured", 
       cex.lab=1.3, cex.main=1.3,
       main="Objective space")
  
  tmp_pit <- sandpit
  freqtab <- as.data.frame(ftable(unlist(pool[,4:(nselect + 3)]))) # check you know what this is?
  tmp_pit[as.numeric(as.character(freqtab$Var1))] <- freqtab$Freq
  tmp_pit[tmp_pit == 0] = NA
  plot(tmp_pit, col=greens(100), xaxt="n", yaxt="n", cex.main=1.3, main="Geographic space")
  mtext("Longitude", 1, 1.5, cex=1.3)
  mtext("Latitude", 2, 1.5, cex=1.3)
  mtext("Initial front", outer=TRUE, side=3, line=2, font=2, cex=1.5)
  dev.off()
  
  starting_front <- pool
  
  for (iter in 1:niters){
    message(iter)
    n_add <- poolsize - nrow(pool)
    
    if (sample_method == "random"){
      additions <- t(sapply(1:n_add, function(i){sample(site_ids, nselect, replace=TRUE)})) %>%
        as.data.frame()
    }
    
    if (sample_method == "neighbours"){
      # under construction
      additions <- t(sapply(0: (n_add-1), function(i){
        neighs <- as.vector(neighbourhood_matrix[unlist(pool[i %% nrow(pool) + 1, 4:(3 + nselect)]),])
        sample(neighs[!is.na(neighs)], size=nselect, replace=FALSE) # not too worried about repeats ... will tend towards designs with less sites which is kinda what I want
      })) %>%
        as.data.frame()
    }
    
    names(additions) <- paste0("site", 1:nselect)
    
    obj <- data.frame(sum_risk = apply(additions, 1, 
                                       function(x){objective_list[[1]](site_ids=unique(as.vector(catchment_matrix[unlist(x),])), 
                                                                       vector_to_be_indexed=potential_vec)}),
                      sum_pop = apply(additions, 1, 
                                      function(x){objective_list[[1]](site_ids=unique(as.vector(catchment_matrix[unlist(x),])), 
                                                                      vector_to_be_indexed=pop_vec)}),
                      addition = rep(TRUE, n_add))
    additions <- cbind(obj, additions)
    pool <- rbind(pool, as.matrix(additions))
    
    png(paste0("~/Desktop/knowlesi/multi_site/output/toy_movie/iter", iter, ".png"), width=2200, height=1400, pointsize=40)
    par(mfrow=c(1,2), mar=c(5.1,4.1,2,2), oma=c(0,0,4.1,1.1))
    
    plot(pool$sum_risk[pool$addition == TRUE], log(pool$sum_pop[pool$addition == TRUE]),
         #main=paste("Step", iter), 
         xlim=box_extent[1:2], ylim=box_extent[3:4],
         xlab="Total potential risk captured", ylab="Total human population captured",
         main="Objective space", cex.main=1.3, cex.lab=1.3)
    points(pool$sum_risk[pool$addition == FALSE], log(pool$sum_pop[pool$addition == FALSE]),
           col=apple, pch=16)
    
    # find pareto front + remove dominated designs
    pool <- psel(pool, high("sum_risk") * high("sum_pop"), top_level=top_level, show_level=FALSE)
    progress[iter,] <- c(sum(pool$addition), nrow(pool), poolsize - nrow(pool))
    
    tmp_pit <- sandpit
    freqtab <- as.data.frame(ftable(unlist(pool[,4:(nselect + 3)])))
    tmp_pit[as.numeric(as.character(freqtab$Var1))] <- freqtab$Freq
    tmp_pit[tmp_pit == 0] = NA
    plot(tmp_pit, col=greens(100), xaxt="n", yaxt="n", main="Geographic space", cex.main=1.3)
    mtext("Longitude", 1, 1.5, cex=1.3)
    mtext("Latitude", 2, 1.5, cex=1.3)
    mtext(paste("Step", iter), outer=TRUE, side=3, line=2, font=2, cex=1.5)
    dev.off()
    sandcastle <- addLayer(sandcastle, tmp_pit)
    
    pool$addition <- rep(FALSE, nrow(pool))
  }
  
  front <- pool[order(pool$sum_risk),]
  front <- starting_front[order(starting_front$sum_risk),]
  
  return(list(pool=pool,
              progress=progress,
              sandcastle=sandcastle))
}