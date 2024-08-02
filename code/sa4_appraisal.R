# script for Victoria-wide appraisal, trimmed down from jev/code/mozzie_surveillance_potential.R

all_mozzies <- read.csv('~/Desktop/jev/from_Freya_local/JEV_secure/data/phase_2/clean/pathogen/all_moz_surveillance.csv') %>%
  rename(host_type=host_species) %>% 
  dplyr::select(
    -c(site,
       n_individual_spp_tests,
    )
  ) %>% 
  drop_na(longitude:virus_detection)
all_mozzies$data_source[all_mozzies$data_source == "VIC"] = "Vic"

potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
hpop <- raster("~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif")

states = st_read("~/Desktop/jev/data/admin/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")
sa4s <- st_read("~/Desktop/jev/data/admin/SA4_2021_AUST_SHP_GDA2020/SA4_2021_AUST_GDA2020.shp")

state_potential_surveillance_compare <- function(state,
                                                 state_shp,
                                                 sa4s,
                                                 potential_national,
                                                 hpop_national,
                                                 past_mozzies,
                                                 hpop_breaks,
                                                 buffer_distance_m=10000,  # for state buffer; should increase if blur radius is increased
                                                 blur_radius_degrees=0.05,  # for site buffer
                                                 quant_prob=0.90,
                                                 out_path="figures/"){
  # here is a big function that should handle everything for us ...
  # I only need it to work for Vic for the chapter but might as well leave the ability
  # to use it elsewhere ...
  state_shadow <- state_shp %>%
    st_buffer(dist = buffer_distance_m) %>%
    st_simplify(dTolerance = 5000)
  
  potential_buffered <- potential_national %>%
    crop(state_shadow) %>%
    raster::mask(state_shadow)
  
  potential_state <- potential_national %>%
    crop(state_shadow) %>%
    raster::mask(state_shp)
  
  hpop_state <- hpop_national %>%
    crop(state_shadow) %>%
    raster::mask(state_shp)
  
  wtmat = focalWeight(potential_buffered, blur_radius_degrees,
                      type="circle")
  wtmat[wtmat != 0] = 1
  
  blurred_potential = focal(potential_buffered, wtmat, sum, na.rm=TRUE)
  blurred_potential_masked = raster::mask(blurred_potential, potential_state)
  
  quant = quantile(blurred_potential_masked, probs=c(quant_prob))
  sel = blurred_potential_masked
  sel[sel < quant] = NA
  
  short_juris = list("New South Wales"="NSW",
                     "Victoria"="Vic",
                     "Queensland"="Qld",
                     "South Australia"="SA",
                     "Western Australia"="WA",
                     "Tasmania"="Tas",
                     "Northern Territory"="NT",
                     "Australian Capital Territory"="ACT")
  thin_juris = c("NT", "Qld", "WA") # the negative space was annoying me ...
  
  state_sa4s <- sa4s[sa4s$STE_NAME21 == names(short_juris)[which(short_juris == state)],]
  state_sa4s <- state_sa4s[!st_is_empty(state_sa4s),] %>%
    st_simplify(dTolerance = 1000)
  
  if (state %in% thin_juris){
    png(paste0(out_path, "surveillance_potential_", state, ".png"),
        height=2700, width=3800, res=300)
  } else {
    png(paste0(out_path, "surveillance_potential_", state, ".png"),
        height=2300, width=4200, res=300)
  }
  
  pn_cols <- colorRampPalette(colors = c("#f7f7f7", "#c23375"))(1000)
  purps = colorRampPalette(brewer.pal(9,"Purples"))(100)
  extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(1000)
  
  par(mfrow=c(1,2), bty="n", mar=c(2.1,1.1,4.1,1.1))
  
  plot(potential_buffered, col=brewer.pal(9, "Greys"), legend=FALSE,
       xaxt="n", yaxt="n", horizontal=TRUE)
  par(new=TRUE, bty="n", mar=c(2.1,1.1,4.1,1.1))
  plot(potential_state, col=pn_cols,  xaxt="n", yaxt="n", horizontal=TRUE,
       legend.args=list("Transmission suitability", side=3, line=1, cex=1.4),
       xaxt="n", yaxt="n",
       axis.args=list(at=c(minValue(potential_state)+0.05, maxValue(potential_state)-0.05),
                      labels=c("Low","High")))
  
  par(new=TRUE, bty="n", mar=c(2.1,1.1,4.1,1.1))
  plot(sel, col=alpha("orange", 0.5), legend=FALSE,  xaxt="n", yaxt="n", horizontal=TRUE)
  par(new=TRUE, bty="n", mar=c(2.1,1.1,4.1,1.1))
  plot(st_geometry(state_shp), add=TRUE, lwd=0.5)
  points(past_mozzies[past_mozzies$data_source == state,
                      c("longitude", "latitude")], 
         pch=0, cex=0.8, col="purple")   
  # legend("topright", fill=alpha("orange", 0.5), cex=1.8,
  #        "Areas of highest\n transmission suitability", bty="n")
  
  plot(sqrt(hpop_state), col=extreme_purps, xaxt="n", yaxt="n", horizontal=TRUE,
       legend.args=list("Distance-weighted human population density", side=3, 
                        line=1, cex=1.4),
       axis.args=list(at=sqrt(hpop_breaks), labels=hpop_breaks))
  plot(sel, col=alpha("orange", 0.5), legend=FALSE, add=TRUE)
  plot(st_geometry(state_sa4s), add=TRUE, lwd=0.5)
  points(past_mozzies[past_mozzies$data_source == state,
                      c("longitude", "latitude")], 
         pch=0, cex=0.8, col="purple")
  
  par(mfrow=c(1,1), new=TRUE, mar=c(2.1,1.1,0,1.1),xpd=NA)
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  legend(0.2,1.05, fill=alpha("orange", 0.5), cex=1.6,
         "Areas of highest \ntransmission suitability", bty="n")
  dev.off()
  
  return(list(potent=potential_state,
              hpop=hpop_state,
              at_risk=sel))
}


comments_col = function(surveil_df, rank_cutoff=12){
  # applies commenting rules to surveil_df (called by surveil_df function)
  par_trap_ranks = data.frame(rank(surveil_df$pc_par, ties.method="average"),
                              rank(surveil_df$pc_traps_state, ties.method="average"))
  par_trap_ranks$diff = par_trap_ranks[,1] - par_trap_ranks[,2]
  
  surveil_df = mutate(surveil_df,
                      comments = case_when(pc_traps_region < 60 & total_traps >= 10 ~ "potential to re-target within SA4",
                                           par_trap_ranks$diff >= rank_cutoff ~ "potential to re-distribute to SA4",
                                           par_trap_ranks$diff <= rank_cutoff*-1 ~ "potential to re-distribute from SA4",
                                           .default = "-"))
  
  return(surveil_df$comments)
}


surveil_df <- function(state,
                       surveil_rasts,
                       sa4s,
                       state_shp,
                       hpop_national,
                       all_traps,
                       rank_cutoff_state){
  # big function to cook up the table ....
  # would like columns for total SA4 population and geographic area ....
  
  hpop_state <- hpop_national %>%
    crop(surveil_rasts$potent) %>%
    mask(surveil_rasts$potent)
  
  surveil_rasts$at_risk[!is.na(surveil_rasts$at_risk)] = 1
  
  state_par = sum(values(hpop_state * surveil_rasts$at_risk), na.rm=TRUE)
  message(paste0("State PAR: ", state_par))
  
  state_sa4s <- sa4s[sa4s$STE_NAME21 == state,]
  state_sa4s <- state_sa4s[!st_is_empty(state_sa4s),]
  state_sa4_ras <- rasterize(state_sa4s, hpop_state, background=NA)
  
  sa4_par <- zonal(hpop_state * surveil_rasts$at_risk, state_sa4_ras, fun = "sum")
  sa4_pop <- zonal(hpop_state, state_sa4_ras, fun = "sum")
  
  # now for the traps ...
  all_traps$at_risk <- raster::extract(surveil_rasts$at_risk, 
                                       all_traps[,c("longitude","latitude")])
  all_traps$sa4 <- raster::extract(state_sa4_ras,
                                   all_traps[,c("longitude","latitude")])
  
  trap_summary <- all_traps[, c("at_risk", "sa4")] %>%
    drop_na(sa4) %>%
    mutate(at_risk=case_when(is.na(at_risk) ~ 0,
                             at_risk == 1 ~ 1)) %>%
    group_by(sa4) %>%
    summarise(traps_at_risk = sum(at_risk), total_traps=n())
  
  short_juris = list("New South Wales"="NSW",
                     "Victoria"="Vic",
                     "Queensland"="Qld",
                     "South Australia"="SA",
                     "Western Australia"="WA",
                     "Tasmania"="Tas",
                     "Northern Territory"="NT",
                     "Australian Capital Territory"="ACT")
  
  out = data.frame(sa4_name=state_sa4s$SA4_NAME21, 
                   par=sa4_par[,2], 
                   total_pop=sa4_pop[,2],
                   traps_at_risk=0,
                   total_traps=0)
  
  out[trap_summary$sa4,
      c("traps_at_risk", "total_traps")] = trap_summary[,c("traps_at_risk", "total_traps")]
  state_traps <- sum(out$total_traps)
  message(paste0("State traps: ", state_traps))
  
  # fix up percentages:
  fix_pc = function(pc){case_when(is.na(pc) ~ "-",
                                  pc == 0 | pc == 100 ~ paste0(pc, "%"),
                                  pc < 1 & pc > 0 ~ "<1%",
                                  pc > 99 & pc < 100 ~ ">99%",
                                  pc >= 1 & pc <= 99 ~ paste0(format(round(pc, digits=0)), "%"))}
  
  # have an "out_numeric" and an "out": "out" is for nice formatting, 
  # "out_numeric" is for applying commenting rules to
  message(paste(names(out), collapse=" "))
  out_numeric = out %>%
    mutate(pc_par = par / state_par * 100,
           pc_traps_state = total_traps / state_traps * 100,
           pc_traps_region = case_when(total_traps == 0 ~ NA,
                                       total_traps != 0 ~ traps_at_risk / total_traps * 100)) #%>%
  
  out_numeric$comments = comments_col(out_numeric, rank_cutoff = rank_cutoff_state)
  
  out = out_numeric %>%
    mutate(#par=format(round(par, digits=0), big.mark = ",", scientific = FALSE),
      par=format(round(par/1000, digits=0)),
      pc_par=fix_pc(pc_par),
      pc_traps_state=fix_pc(pc_traps_state),
      pc_traps_region=fix_pc(pc_traps_region),
      total_traps=format(round(total_traps, digits=0)),
      traps_at_risk=format(round(traps_at_risk, digits=0)))
  
  # adding total pop in here?
  out = out[,c("sa4_name","total_pop","par","total_traps", "pc_par",
               "pc_traps_state", "traps_at_risk","pc_traps_region", "comments")]
  names(out) = c("SA4",
                 "Total pop.",
                 "PAR", 
                 "No' catches", 
                 paste0("% ", short_juris[[state]], " PAR"), 
                 paste0("% ", short_juris[[state]], " catches"),
                 "No' catches in risk", 
                 "% catches in risk", 
                 "Comments")
  
  return(list(out_numeric=out_numeric,
              out=out))
}




# make up some results for vic .....


vic_shp <- states %>%
  dplyr::filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)
vic_surveil <- state_potential_surveillance_compare("Vic", vic_shp, sa4s,
                                                    potential, 
                                                    hpop, all_mozzies,
                                                    hpop_breaks=c(25,100, 250, 500))

# play around with rank_cutoff here?
vic_surveil_df = surveil_df("Victoria",
                            vic_surveil,
                            sa4s,
                            vic_shp,
                            hpop,
                            all_mozzies,
                            9)

library(xtable)
print(xtable::xtable(vic_surveil_df$out),
      include.rownames=FALSE)



