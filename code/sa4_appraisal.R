# script for Victoria-wide appraisal, trimmed down from jev/code/mozzie_surveillance_potential.R

potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
hpop <- raster("~/Desktop/jev/from_Freya_local/JEV/output/hpop_blur_aus_0.2_res_0.008.tif")

states = st_read("~/Desktop/jev/data/admin/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")
sa4s <- st_read("~/Desktop/jev/data/admin/SA4_2021_AUST_SHP_GDA2020/SA4_2021_AUST_GDA2020.shp")
sa4s$SA4_NAME21[sa4s$STE_NAME21 == "Victoria"]

vic_shp <- states %>%
  dplyr::filter(STE_NAME21 == "Victoria") %>%
  st_simplify(dTolerance = 1000)

vic_sa4s <- sa4s[sa4s$STE_NAME21 == "Victoria",] %>%
  st_simplify(dTolerance = 1000)

wa_shp <- states %>%
  dplyr::filter(STE_NAME21 == "Western Australia") %>%
  st_simplify(dTolerance = 1000)

wa_sa4s <- sa4s[sa4s$STE_NAME21 == "Western Australia",] %>%
  st_simplify(dTolerance = 1000)

lab_deets_vic <- do.call(rbind, st_geometry(st_centroid(vic_sa4s))) %>%
  as.data.frame() %>%
  setNames(c("centlon", "centlat")) %>%
  drop_na() %>%
  mutate(name = vic_sa4s$SA4_NAME21[1:17] %>%
           gsub(pattern="Melbourne - North West", replacement="North West\n(Melbourne)") %>%
           gsub(pattern="Melbourne - ", replacement="") %>%
           gsub(pattern=" - ", replacement="-") %>%
           gsub(pattern=" and ", replacement="\nand "),
         lablon = c(142.5, 144.2, 143.2, 148.5, 149.1, 144.9, 148.7, 147, 147, 
                    145, 148.5, 148, 144.2, 146.3, 143, 146, 141),
         lablat = c(-38.7, -35, -39, -36, -36.5, -39.2, -38.5, -39.1, -35.5, 
                    -34.5, -38.2, -38.8, -38.7, -39.4, -34, -35, -39),
         pos = c(1,3,1,3,3,1,4,4,3,3,4,4,1,4,3,3,1))

lab_deets_wa <- do.call(rbind, st_geometry(st_centroid(wa_sa4s))) %>%
  as.data.frame() %>%
  setNames(c("centlon", "centlat")) %>%
  drop_na() %>%
  mutate(name = wa_sa4s$SA4_NAME21[1:10] %>%
           gsub(pattern="Perth - ", replacement="") %>%
           gsub(pattern="Western Australia - ", replacement=""),
         lablon = c(rep(113,7),126,115, 114),
         lablat = c(-35, -34.5, -32.5, -32,-31.5,-33.5,-33,-34,-18,-19),
         pos = c(2,2,2,2,2,2,2,4,2,2))


# this is the mozzie data I use in JEV 1
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
  filter(as.Date(date) > as.Date("2022-07-01"))

wa_mozzies <- read.csv('~/Desktop/jev/from_Freya_local/JEV_secure/data/phase_2/clean/pathogen/all_moz_surveillance.csv') %>%
  rename(host_type=host_species) %>% 
  dplyr::select(
    -c(site,
       n_individual_spp_tests,
    )
  ) %>% 
  drop_na(longitude:virus_detection) %>%
  subset(data_source == "WA") %>%
  filter(as.Date(date) > as.Date("2022-07-01"))

# number of unique locations
nrow(unique(vic_mozzies[,c("longitude","latitude")]))
nrow(unique(wa_mozzies[,c("longitude","latitude")]))

# date range, number of trapping events
range(vic_mozzies$date)
range(wa_mozzies$date)
nrow(vic_mozzies)
nrow(wa_mozzies)

# number of unique pixels (fine case)
id_ras <- potential
values(id_ras) <- 1:ncell(potential)
vic_mozzies %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise() %>%
  nrow()
wa_mozzies %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise() %>%
  nrow()

# number of unique pixels (coarse case)
id_ras <- potential %>%
  aggregate(AGG_FACTOR)
values(id_ras) <- 1:ncell(id_ras)
vic_mozzies %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise() %>%
  nrow()
wa_mozzies %>%
  mutate(pix = raster::extract(id_ras, cbind(longitude, latitude))) %>%
  group_by(pix) %>%
  summarise() %>%
  nrow()

# function straight out of the report: 
# give me a map of objectives with existing surveillance overlaid
state_potential_surveillance_compare <- function(state,
                                                 state_shp,
                                                 sa4s,
                                                 potential_national,
                                                 hpop_national,
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
  
  hpop_buffered <- hpop_national %>%
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
  # state_sa4s <- state_sa4s[!st_is_empty(state_sa4s),] %>%
  #   st_simplify(dTolerance = 1000)
  
  # if (state %in% thin_juris){
  #   png(paste0(out_path, "surveillance_potential_", state, ".png"),
  #       height=2700, width=3800, res=300)
  # } else {
  #   png(paste0(out_path, "surveillance_potential_", state, ".png"),
  #       height=2300, width=4200, res=300)
  # }
  
  # par(mfrow=c(1,2), bty="n", mar=c(2.1,1.1,4.1,1.1))
  
  # par(mfrow=c(1,1), new=TRUE, mar=c(2.1,1.1,0,1.1),xpd=NA)
  # plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  # legend(0.78,1.05, fill=alpha("orange", 0.5), cex=1.6,
  #        "Areas of highest \ntransmission suitability", bty="n")
  # subfigure_label(par()$usr, 0, 0.9, "(a)", 1.5)
  # subfigure_label(par()$usr, 0.5, 0.9, "(b)", 1.5)
  # dev.off()
  
  return(list(potent=potential_state,
              potential_buffered=potential_buffered,
              hpop_buffered=hpop_buffered,
              hpop=hpop_state,
              at_risk=sel))
}

vic_surveil <- state_potential_surveillance_compare("Vic", 
                                                    vic_shp, 
                                                    vic_sa4s,
                                                    potential,
                                                    hpop,
                                                    blur_radius_degrees = 0.12) # to line up with coarse scale

wa_surveil <- state_potential_surveillance_compare("WA", 
                                                   wa_shp, 
                                                   wa_sa4s,
                                                   potential, 
                                                   hpop,
                                                   blur_radius_degrees = 0.12)

pn_cols <- colorRampPalette(colors = c("#f7f7f7", "#c23375"))(1000)
purps = colorRampPalette(brewer.pal(9,"Purples"))(100)
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(1000)

# have popped the plotting outside the function because I was getting fed! up!
{png("figures/surveillance_potential_vic_wa.png",
           height=4800, width=4000, res=300)
  mar=c(2.1,1.1,1.1,1.1)
  oma=c(16,2,0,0)
par(mfrow=c(2,2), mar=mar, oma=oma, bty="n")

plot(vic_surveil$potential_buffered, 
     col=colorRampPalette(brewer.pal(9, "Greys"))(100), 
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(1,1))
plot(vic_surveil$potent, col=pn_cols,
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(1,1))
plot(vic_surveil$at_risk, col=alpha("orange", 0.5),
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(1,1))
plot(st_geometry(vic_shp), add=TRUE, lwd=0.5)
points(vic_mozzies[,c("longitude","latitude")], pch=0, cex=0.8, col="purple")
plot(vic_surveil$potent, legend.only=TRUE,
     legend.args=list("Transmission suitability", side=1, line=3, cex=1.4),
     cex.axis=1.2, legend.width=1.5, horizontal=TRUE,
     col=pn_cols)

par(xpd=NA)
xpos=149.75
ypos=-40.5
ycent=-36.5
axis(1, at=c(xpos, xpos + 100/(111.320*cos(ycent/180))), pos=ypos, labels=c("",""))
text(xpos + 100/(111.320*cos(ycent/180))/2, ypos-0.5, "100 km", cex=1.5)

plot(sqrt(vic_surveil$hpop_buffered), 
     col=colorRampPalette(brewer.pal(9, "Greys"))(100),
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(1,2))
hpop_breaks=c(25,100, 250, 500)
plot(sqrt(vic_surveil$hpop), col=purps(100),  
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
plot(sqrt(vic_surveil$hpop), legend.only=TRUE,
     legend.args=list("Distance-weighted human population density", side=1, line=3, cex=1.4),
     axis.args=list(at=sqrt(hpop_breaks), labels=hpop_breaks, cex=1.2),
     cex.axis=1.2, legend.width=1.5, horizontal=TRUE,
     col=purps(100))
par(new=TRUE, mar=mar, oma=oma, mfg=c(1,2))
plot(vic_surveil$at_risk, col=alpha("orange", 0.5),
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(1,2))
plot(st_geometry(vic_sa4s), add=TRUE, lwd=0.5)
points(vic_mozzies[,c("longitude","latitude")], pch=0, cex=0.8, col="purple")

par(xpd=NA)
text(lab_deets_vic$lablon, lab_deets_vic$lablat, labels=lab_deets_vic$name, 
     offset=0.5, pos=lab_deets_vic$pos,
     cex=1.1)
for (i in 1:nrow(lab_deets_vic)){
  lines(lab_deets_vic[i, c("lablon", "centlon")], 
        lab_deets_vic[i, c("lablat", "centlat")], lwd=0.7)
}


# reset margins for WA (she's a little longer than wider)
mar=c(4.1,1.1,0.1,1.1)
oma=c(6,2,0,0)
par(mfrow=c(2,2), mar=mar, oma=oma, new=TRUE, mfg=c(2,1), bty="n", xpd=NA)

plot(wa_surveil$potential_buffered, 
     col=colorRampPalette(brewer.pal(9, "Greys"))(100), 
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(2,1))
plot(wa_surveil$potent, col=pn_cols,
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(2,1))
plot(wa_surveil$at_risk, col=alpha("orange", 0.5),
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(2,1), xpd=NA)
plot(st_geometry(wa_shp), add=TRUE, lwd=0.5)
points(wa_mozzies[,c("longitude","latitude")], pch=0, cex=0.8, col="purple")

par(xpd=NA)
xpos=131
ypos=-36.7
ycent=-25
axis(1, at=c(xpos, xpos + 100/(111.320*cos(ycent/180))), pos=ypos, labels=c("",""))
text(xpos + 100/(111.320*cos(ycent/180))/2, ypos-1, "100 km", cex=1.5)


par(mfrow=c(2,2), mar=mar, oma=oma, mfg=c(2,2), new=TRUE)
plot(sqrt(wa_surveil$hpop_buffered), 
     col=colorRampPalette(brewer.pal(9, "Greys"))(100),
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(2,2))
hpop_breaks=c(1,25,100,300)
plot(sqrt(wa_surveil$hpop), col=purps(100),  
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(2,2))
plot(wa_surveil$at_risk, col=alpha("orange", 0.5),
     legend=FALSE, xaxt="n", yaxt="n", horizontal=TRUE, legend.mar=0)
par(new=TRUE, mar=mar, oma=oma, mfg=c(2,2), xpd=NA)
plot(st_geometry(wa_sa4s), add=TRUE, lwd=0.5)
points(wa_mozzies[,c("longitude","latitude")], pch=0, cex=0.8, col="purple")

par(xpd=NA)
text(lab_deets_wa$lablon, lab_deets_wa$lablat, labels=lab_deets_wa$name, 
     offset=0.5, pos=lab_deets_wa$pos,
     cex=1.1)
for (i in 1:nrow(lab_deets_wa)){
  lines(lab_deets_wa[i, c("lablon", "centlon")], 
        lab_deets_wa[i, c("lablat", "centlat")], lwd=0.7)
}

# can I fenangle the WA legends a little lower?
mar=c(2.1,1.1,0.1,1.1)
oma=c(4,2,0,0)
par(mfrow=c(2,2), mar=mar, oma=oma, new=TRUE, mfg=c(2,1))
plot(wa_surveil$potent, legend.only=TRUE,
     legend.args=list("Transmission suitability", side=1, line=3, cex=1.4),
     cex.axis=1.2, legend.width=1.5, horizontal=TRUE,
     col=pn_cols)
par(mfrow=c(2,2), mar=mar, oma=oma, new=TRUE, mfg=c(2,2))
plot(sqrt(wa_surveil$hpop), legend.only=TRUE,
     legend.args=list("Distance-weighted human population density", side=1, line=3, cex=1.4),
     axis.args=list(at=sqrt(hpop_breaks), labels=hpop_breaks, cex=1.2),
     cex.axis=1.2, legend.width=1.5, horizontal=TRUE,
     col=purps(100))



# Finishing touches ...
par(mfrow=c(1,1), new=TRUE, mar=c(2.1,1.1,0,1.1), xpd=NA)
plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
legend(0.78,1.05, fill=alpha("orange", 0.5), cex=1.5,
       "Areas of highest \ntransmission suitability", bty="n")
subfigure_label(par()$usr, 0, 0.95, "(a)", 1.5)
subfigure_label(par()$usr, 0.5, 0.95, "(b)", 1.5)
subfigure_label(par()$usr, 0, 0.5, "(c)", 1.5)
subfigure_label(par()$usr, 0.5, 0.5, "(d)", 1.5)

dev.off()}

# which are hitting ?
# plot(sqrt(wa_surveil$hpop_buffered), xlim=c(113,120), ylim=c(-35,-30),
#      col=colorRampPalette(brewer.pal(9, "Greys"))(100),
#      legend=FALSE, horizontal=TRUE, legend.mar=0)
# hpop_breaks=c(1,25,100,300)
# plot(sqrt(wa_surveil$hpop), col=purps(100), add=TRUE,
#      legend=FALSE)
# plot(wa_surveil$at_risk, col=alpha("orange", 0.5), add=TRUE,
#      legend=FALSE, horizontal=TRUE, legend.mar=0)
# plot(st_geometry(wa_sa4s), add=TRUE, lwd=0.5)
# tmp <- raster::extract(wa_surveil$at_risk, wa_mozzies[,c("longitude","latitude")])
# points(wa_mozzies[!is.na(tmp),c("longitude","latitude")], pch=4, cex=0.8, col="purple")
# points(wa_mozzies[is.na(tmp),c("longitude","latitude")], pch=4, cex=0.8, col="red")

# I don't reckon there's enough room for an inset figure here
{png("figures/surveillance_potential_swwa.png",
    height=1950,
    width=2200,
    pointsize=45)
par(mar=c(5.1,4.1,4.1,10.1))
plot(wa_surveil$potent, xlim=c(114,120), ylim=c(-35.1,-30), col=pn_cols,
     bty="n", xlab="Longitude", ylab="Latitude", legend.mar=0, legend=FALSE,
     main="Transmission suitability: south-west WA", cex.lab=1.2, cex.main=1.2)
     #legend.args=list(text ="Transmission suitability", side=4))
plot(wa_surveil$at_risk, col=alpha("orange", 0.5), add=TRUE, legend=FALSE)
par(xpd=FALSE)
plot(st_geometry(wa_sa4s), lwd=2, add=TRUE)
points(wa_mozzies[,c("longitude","latitude")], pch=0, col="purple", lwd=3)
par(mar=c(5.1,4.1,4.1,2.1), new=TRUE)
plot(wa_surveil$potent, xlim=c(114,120), ylim=c(-35.1,-30), col=pn_cols,
     legend.only=TRUE, legend.args=list(text ="Transmission suitability", side=4, line=-3, cex=1.2))
par(mfrow=c(1,1), new=TRUE, mar=c(2.1,1.1,0,1.1), xpd=NA, bty="n")
plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
legend(0.8, 0.2, "Areas of highest\ntransmission\nsuitability", 
       fill=alpha("orange", 0.5), bty="n", cex=1.1)
dev.off()}

# it's probably time to pop the plotting out so that I can do these two in one figure

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
                   sa4_area=state_sa4s$AREASQKM21,
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
      sa4_area=format(round(sa4_area, digits=0), big.mark = ","),
      total_pop=format(round(total_pop/1000, digits=0)),
      par=format(round(par/1000, digits=0)),
      pc_par=fix_pc(pc_par),
      pc_traps_state=fix_pc(pc_traps_state),
      pc_traps_region=fix_pc(pc_traps_region),
      total_traps=format(round(total_traps, digits=0)),
      traps_at_risk=format(round(traps_at_risk, digits=0)))
  
  # adding total pop in here?
  out = out[,c("sa4_name","sa4_area","total_pop","par","total_traps", "pc_par",
               "pc_traps_state", "traps_at_risk","pc_traps_region", "comments")]
  names(out) = c("SA4",
                 "Area (kmsq)",
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



# play around with rank_cutoff here?
vic_surveil_df = surveil_df("Victoria",
                            vic_surveil,
                            vic_sa4s,
                            vic_shp,
                            hpop,
                            vic_mozzies,
                            9)

wa_sa4s$SA4_NAME21 <- gsub("Western Australia - ", "", wa_sa4s$SA4_NAME21)
wa_surveil_df = surveil_df("Western Australia",
                           wa_surveil,
                           wa_sa4s,
                           wa_shp,
                           hpop,
                           wa_mozzies,
                           6)


library(xtable)

add.to.row <- list(pos=list(-1),
                   command=c("\\multicolumn{2}{c}{Group 1 and 2} & \\multicolumn{1}{c}{Group 3} \\\\ \\cmidrule(lr){1-2} \\cmidrule(lr){3-3}"))
print(xtable::xtable(vic_surveil_df$out),
      include.rownames=FALSE)

print(xtable::xtable(wa_surveil_df$out),
      include.rownames=FALSE)



