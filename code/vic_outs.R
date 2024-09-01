# process vic outputs and run up vic_sensitivity.png and vic_mapped.png
# need to run what's in main.R
source("code/vic_setup.R")

vic_mozzies <- read.csv('~/Desktop/jev/from_Freya_local/JEV_secure/data/phase_2/clean/pathogen/all_moz_surveillance.csv') %>%
  rename(host_type=host_species) %>% 
  dplyr::select(
    -c(site,
       n_individual_spp_tests,
    )
  ) %>% 
  drop_na(longitude:virus_detection) %>%
  mutate(data_source = case_when(data_source == "VIC" ~ "Vic",
                                 TRUE ~ data_source)) %>%
  subset(data_source == "Vic") %>%
  filter(as.Date(date) > as.Date("2022-07-01")) %>%
  group_by(longitude, latitude) %>%
  summarise(ntimes = n()) %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise()

objective_func = function(x, catch_mem, vec){
  sum(vec[unique(as.vector(catch_mem[unlist(x),]))], na.rm=TRUE)
}

# not sure if I believe how good this is under risk obj ... although there are more sites here
existing_potent <- objective_func(vic_mozzies$pix, 
                                  catch_membership_mat, 
                                  vic_objective$buffer_potent)
existing_hpop <- objective_func(vic_mozzies$pix,
                                catch_membership_mat,
                                vic_objective$buffer_hpop)

nrow(vic_mozzies)


progress1000 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10.csv")
progress5000 <- read.csv("output/vic/vic_auc_pool5000_iters100_runs10.csv")
progress10000 <- read.csv("output/vic/vic_auc_pool10000_iters100_runs10.csv")
progress50000 <- read.csv("output/vic/vic_auc_pool50000_iters100_runs10.csv")

progress_neigh1 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10.csv")
progress_neigh2 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_neigh2.csv")
progress_neigh3 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_neigh3.csv")
progress_neigh4 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_neigh4.csv")
progress_neigh5 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_neigh5.csv")

progress_pareto1 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10.csv")
progress_pareto2 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_pareto2.csv")
progress_pareto3 <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_pareto3.csv")

progress_uneducated <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10.csv")
progress_educated <- read.csv("output/vic/vic_auc_pool1000_iters100_runs10_greedystart.csv")

progress_apple <- read.csv("output/vic/vic_auc_pool50000_iters100_runs10_neigh3.csv")
progress_pear <- read.csv("output/vic/vic_auc_pool50000_iters100_runs10_neigh3_greedystart.csv")

#progress_educated <- read.csv("output/old_rasters/vic_auc_pool1000_iters100_runs10_greedystart.csv")
#load("output/old_rasters/diagnostics_vic_greedy.rds")

load("output/vic/diagnostics_vic_pool.rds")
load("output/vic/diagnostics_vic_neigh.rds")
load("output/vic/diagnostics_vic_pareto.rds")
load("output/vic/diagnostics_vic_greedy.rds")
load("output/vic/diagnostics_vic_50000_neigh3_greedy.rds")
load("output/vic/diagnostics_vic_50000_neigh3.rds")

# incorporate the rest of the bits in once I've got them and set common ylim
{png("figures/vic_sensitivity.png",
     width=2400,
     height=2000,
     pointsize = 40)
  
  par(mfrow=c(2,2), mar=c(2.1,2.1,2.1,2.1), oma=c(3,4,1,0))
  ylim=c(3397865,6691487)
  
  pool_lim <- auc_agg_fig(list(progress1000,
                               progress5000,
                                progress10000,
                                progress50000),
                          legend_labs=c("1,000", "5,000", "10,000", "50,000"),
                          legend_title="Pool size",
                          pal=c(iddu(4)[2:4], brat),
                          ylim=ylim)
  
  neigh_lim <- auc_agg_fig(list(progress_neigh1,
                                progress_neigh2,
                                progress_neigh3,
                                progress_neigh4,
                                progress_neigh5),
                           #progress_neigh10),
                           legend_labs=c("1-neighbours", "2-neighbours", 
                                         "3-neighbours", "4-neighbours","5-neighbours"),# "10-neighbours"),
                           legend_title="Neighbourhood size",
                           pal=c(iddu(4)[2:4], brat, iddu(4)[1]),
                           ylim=ylim)
  
  pareto_lim <- auc_agg_fig(list(progress_pareto1,
                                 progress_pareto2,
                                 progress_pareto3),
                            legend_labs=c("Rank 1", "Rank 2", "Rank 3"),
                            legend_title="Goldberg ranking cutoff",
                            pal=iddu(4)[2:4],
                            ylim=ylim)
  
  loaded_lim <- auc_agg_fig(list(progress_uneducated,
                                 progress_educated),
                            legend_labs=c("Random", "Greedy"),
                            legend_title="Starting pool",
                            pal=iddu(4)[2:4],
                            ylim=ylim)
  
  mtext("Area under estimated Pareto front", 2, outer=TRUE, line=2, cex=1.2)
  mtext("Iteration", 1, outer=TRUE, line=1, cex=1.2)
  
  par(new=TRUE, mfrow=c(1,1), xpd=NA)
  empty_plot_for_legend()
  subfigure_label(par()$usr, -0.05, 1.06, "(a)")
  subfigure_label(par()$usr, 0.515, 1.06, "(b)")
  subfigure_label(par()$usr, -0.05, 0.48, "(c)")
  subfigure_label(par()$usr, 0.515, 0.48, "(d)")
  dev.off()}

range(pool_lim, neigh_lim, pareto_lim, loaded_lim)
#################################################################################
# vic map figure
message("Make sure this is our best guess: ")
final_frontsdf <- final_fronts_apple %>%
  append(final_fronts_pear) %>%
  append(final_fronts_neigh3) %>% 
  append(final_frontsgreedy) %>%
  append(final_fronts50000) %>%
  rbindlist() %>%
  as.data.frame()
agg_pareto <- psel(final_frontsdf, high("sum_pop")*high("sum_risk")) %>%
  arrange(sum_pop)

# a little concerned the same design is in here a bunch of times?
all_sites <- agg_pareto %>%
  dplyr::select(grep("site", names(agg_pareto))) %>%
  unique() %>%
  unlist() %>%
  ftable() %>%
  as.data.frame()

give_me_design_and_catch <- function(ras, catch_mat, df, row_ind){
  # for the _mapped figures
  values(ras) <- NA
  catch <- ras
  values(ras)[unlist(df[row_ind,grep("site", names(df))])] <- 1
  values(catch)[unique(as.vector(catch_mat[unlist(df[row_ind,grep("site", names(df))]),]))] <- 1
  return(list(design=ras,
              catch=catch))
}

pot_razzes <- give_me_design_and_catch(vic_objective$potent,
                                       catch_membership_mat,
                                       agg_pareto,
                                       1)

pop_razzes <- give_me_design_and_catch(vic_objective$potent,
                                       catch_membership_mat,
                                       agg_pareto,
                                       nrow(agg_pareto))

mid <- nrow(agg_pareto)/2
mid_razzes <- give_me_design_and_catch(vic_objective$potent,
                                       catch_membership_mat,
                                       agg_pareto,
                                       mid)

agg_map <- vic_objective$potent
values(agg_map)[!is.na(values(agg_map))] <- 0
values(agg_map)[as.numeric(paste(all_sites$.))] <- all_sites$Freq
values(agg_map)[values(agg_map) == 0] <- NA

# existing surveillance
actual_map <- vic_objective$potent
values(actual_map) <- NA
values(actual_map)[vic_mozzies$pix] <- 1

actual_catch <- vic_objective$potent
values(actual_catch) <- NA
values(actual_catch)[unique(as.vector(catch_membership_mat[vic_mozzies$pix,]))] <- 1

{png("figures/vic_mapped.png",
     width=2500,
     height=2000,
     pointsize=35)
  par(mfrow=c(1,1), oma=c(18,0,2,24), mar=c(4,4,2,0))
  
  plot(final_frontsdf$sum_pop, final_frontsdf$sum_risk, 
       xlab="Sum(Human Population Density)",
       ylab="Sum(Potential Risk)",
       cex.lab=1.2, cex.axis=1.2,
       xlim=range(existing_hpop, final_frontsdf$sum_pop))
  lines(rep(max(agg_pareto$sum_pop), 2), c(0, agg_pareto$sum_risk[nrow(agg_pareto)]), col="grey", lty=2, lwd=2)
  lines(c(0, agg_pareto$sum_pop[1]), rep(max(agg_pareto$sum_risk), 2), col="grey", lty=2, lwd=2)
  points(agg_pareto$sum_pop, agg_pareto$sum_risk, col=brat, pch=16)
  lines(agg_pareto$sum_pop, agg_pareto$sum_risk, col=brat, lwd=2)
  #text(30000, 70, paste("AUF:\n", format(round(pareto_progress_auc(list(agg_pareto)), digits=0), big.mark=",")),
  #     col=alpha(brat, 0.6), cex=2.5, font=2)
  text(existing_hpop, 135, "Existing\n surveillance", col=iddu(2)[2], cex=1.2)
  #arrows(existing_hpop, 93, existing_hpop, existing_potent-1, col=iddu(2)[2], lwd=2)
  points(agg_pareto[c(1,mid,nrow(agg_pareto)), c("sum_pop", "sum_risk")], 
         col=c(berry, "orange", purp), pch=16, cex=1.5)
  points(agg_pareto[c(1,mid,nrow(agg_pareto)), c("sum_pop", "sum_risk")], 
         col=c(berry, "orange", purp), cex=3, lwd=2)
  points(existing_hpop, existing_potent, col=iddu(2)[2], pch=16, cex=1.5) # make this a little easier to see?
  points(existing_hpop, existing_potent, col=iddu(2)[2], cex=3, lwd=2)
  
  # make it a 3*3 !
  par(mfrow=c(3,3), oma=c(0,2,0,0), mar=c(0,0,0,0), mfg=c(1,3), bty="n", new=TRUE)
  plot(agg_map, col=greens(100), axes=FALSE, bty="n", 
       legend=FALSE)
  plot(agg_map, col=greens(100), legend.only=TRUE, 
       horizontal=TRUE,
       #smallplot= c(-0.05, -0.05 + 0.015, 0.35, 0.7),
       legend.args=list(text="Frequency selected", side=1, cex=1.1, line=2.5)) # could put this in the top right ..
  plot(st_geometry(vic_shp), add=TRUE)
  
  par(mfrow=c(3,3), oma=c(0,2,0,0), mar=c(0,0,0,0), mfg=c(2,3), bty="n", new=TRUE)
  # at least one of these look like they're over the border?
  plot(actual_catch, col=alpha(iddu(4)[2], 0.2), axes=FALSE, bty="n", legend=FALSE)
  plot(actual_map, col=iddu(4)[2], bty="n", legend=FALSE, add=TRUE)
  plot(st_geometry(vic_shp), add=TRUE)
  
  par(mfrow=c(3,3), oma=c(0,2,0,0), mar=c(0,0,0,0), new=TRUE, mfg=c(3,1))
  plot(pot_razzes$catch, col=alpha(berry, 0.2), axes=FALSE, bty="n", legend=FALSE) # weird that there are overlapping catchments in here ...
  plot(pot_razzes$design, add=TRUE, col=berry, legend=FALSE)
  # add map of catchment?
  plot(st_geometry(vic_shp), add=TRUE)
  
  plot(mid_razzes$catch, col=alpha("orange", 0.2), axes=FALSE, bty="n", legend=FALSE)
  plot(mid_razzes$design, col="orange", add=TRUE, legend=FALSE)
  plot(st_geometry(vic_shp), add=TRUE, legend=FALSE)
  
  plot(pop_razzes$catch, col=alpha(purp, 0.2), axes=FALSE, bty="n", legend=FALSE)
  plot(pop_razzes$design, col=purp, add=TRUE, legend=FALSE) # or plot points as in the toy fig ...
  plot(st_geometry(vic_shp), add=TRUE)
  
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(2,2,2,2), new=TRUE, bty="o", xpd=NA)
  empty_plot_for_legend()
  subfigure_label(par()$usr, 0,1,"(a)", 1.2)
  subfigure_label(par()$usr, 0.68,1,"(b)", 1.2)
  subfigure_label(par()$usr, 0.68,0.63,"(c)", 1.2)
  subfigure_label(par()$usr, -0.02,0.28,"(d)", 1.2)
  subfigure_label(par()$usr, 0.33,0.28,"(e)", 1.2)
  subfigure_label(par()$usr, 0.68,0.28,"(f)", 1.2)
  dev.off()}




