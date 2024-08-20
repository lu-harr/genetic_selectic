# methods figure: show objective surfaces
{png("figures/toy_problem_obj.png",
     height=1000,
     width=2200,
     pointsize=40)
  par(mfrow=c(1,2), mar=c(0.1,1.1,4.1,6.1))
  plot(guelphia_potential, col=pinks(100), main="Potential risk of \nJEV transmission", 
       axes=FALSE, xlab="", ylab="", legend=FALSE, cex.main=1.4)
  plot(guelphia_potential, col=pinks(100), legend.only=TRUE,
       legend.args=list(text="Potential risk", side=4, line=-2, cex=1.2),
       legend.width=1.2)
  # actually may not need to log this ...
  plot(guelphia_hpop, col=purps(100), main="Human population density",
       axes=FALSE, xlab="", ylab="", legend=FALSE, cex.main=1.4) # legend will need a title :)
  plot(guelphia_hpop, col=purps(100), legend.only=TRUE,
       legend.args=list(text="Human population density", side=4, line=-2, cex=1.2),
       legend.width=1.2) # units?
  
  par(xpd=NA, new=TRUE, mfrow=c(1,1))
  empty_plot_for_legend()
  subfigure_label(par()$usr, 0, 1.12, "(a)", 1.2)
  subfigure_label(par()$usr, 0.58, 1.12, "(b)", 1.2)
  dev.off()}

############################################################################
# demo ``neighbourhood'' concept for methods: 
# pink transmission suitability bg
# with orange catch overlaid

pn_cols <- colorRampPalette(colors = c("#f7f7f7", "#c23375"))(1000)
ext = c(147.3,147.7,-35.2,-34.8)

potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
pot_eg = crop(potential, extent(ext))
maskmat = pot_eg
values(maskmat) = NA
maskmat = as.matrix(maskmat)
allx = unique(rasterToPoints(pot_eg)[,1])
ally = unique(rasterToPoints(pot_eg)[,2])
mid = c(allx[length(allx)/2], ally[length(ally)/2])

wtmat = focalWeight(pot_eg, 0.05,
                    type="circle")
wtmat[wtmat != 0] = 1

maskmat[((dim(maskmat)[1] - dim(wtmat)[1] + 1)/2):((dim(maskmat)[1] - dim(wtmat)[1] + 1)/2 + dim(wtmat)[1] - 1),
        ((dim(maskmat)[2] - dim(wtmat)[2] + 1)/2):((dim(maskmat)[2] - dim(wtmat)[2] + 1)/2 + dim(wtmat)[2] - 1)] = wtmat
maskmat[maskmat == 0] = NA

maskras = pot_eg
values(maskras) = maskmat

{png("figures/neighbourhood_eg.png",
     height=1800,
     width=2200, pointsize=50)
  par(bty="n", mar=c(4.1,0.1,4.1,4.1))
  plot(pot_eg, col=pn_cols, xaxt="n", yaxt="n",
       axis.args=list(at=c(minValue(pot_eg), maxValue(pot_eg)), 
                      labels=c("Low","High"), cex.axis=1.4),
       legend.args=list("JEV transmission suitability", side=2, line=2, cex=1.4),
       main="Neighbourhood around a point", cex.main=1.4,
       legend.mar=12)
  axis(1, at=c(147.45, 147.55), labels=c("",""), lwd=3)
  mtext("0.1 degrees", 1, 1, cex=1.3)
  plot(maskras, col=alpha("orange", 0.5), add=TRUE, legend=FALSE)
  points(mid[1], mid[2], pch=4, lwd=5, cex=1.5)
  
  par(mfrow=c(1,1), new=TRUE, mar=c(0,0,0,0))
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  legend(0.65,0.9, fill=alpha("orange", 0.5), "Pixels in \nneighbourhood", 
         cex=1.4, bty="n")
  dev.off()}

{png("figures/catchment_eg.png",
     height=1800,
     width=2200, pointsize=50)
  par(bty="n", mar=c(4.1,0.1,4.1,4.1))
  plot(pot_eg, col=pn_cols, xaxt="n", yaxt="n",
       axis.args=list(at=c(minValue(pot_eg), maxValue(pot_eg)), 
                      labels=c("Low","High"), cex.axis=1.4),
       legend.args=list("JEV transmission suitability", side=2, line=2, cex=1.4),
       main="Catchment around a point", cex.main=1.4,
       legend.mar=12)
  axis(1, at=c(147.45, 147.55), labels=c("",""), lwd=3)
  mtext("0.1 degrees", 1, 1, cex=1.3)
  plot(maskras, col=alpha("orange", 0.5), add=TRUE, legend=FALSE)
  points(mid[1], mid[2], pch=4, lwd=5, cex=1.5)
  
  par(mfrow=c(1,1), new=TRUE, mar=c(0,0,0,0))
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  legend(0.65,0.9, fill=alpha("orange", 0.5), "Pixels in \ncatchment", 
         cex=1.4, bty="n")
  dev.off()}


ext = c(147.35,147.65,-35.15,-34.85)
potential <- raster('~/Desktop/jev/from_Freya_local/JEV/output/continuous suit vectors and avian.tif')
pot_eg = crop(potential, extent(ext))
maskmat = pot_eg
values(maskmat) = NA
maskmat = as.matrix(maskmat)
allxy = rasterToPoints(pot_eg)
allx = unique(allxy[,1])
ally = unique(allxy[,2])
mid = c(allx[length(allx)/2], ally[length(ally)/2])

wtmat = focalWeight(pot_eg, 0.05,
                    type="circle")
wtmat[wtmat != 0] = 1

maskmat[((dim(maskmat)[1] - dim(wtmat)[1] + 1)/2):((dim(maskmat)[1] - dim(wtmat)[1] + 1)/2 + dim(wtmat)[1] - 1),
        ((dim(maskmat)[2] - dim(wtmat)[2] + 1)/2):((dim(maskmat)[2] - dim(wtmat)[2] + 1)/2 + dim(wtmat)[2] - 1)] = wtmat
maskmat[maskmat == 0] = NA

maskras = pot_eg
values(maskras) = maskmat

coarse_catch <- guelphia_potential
values(coarse_catch) <- NA
values(coarse_catch)[c(34:36, 44:46, 54:56)] <- 1


{png("figures/catchment_eg.png",
     height=1700,
     width=2200, pointsize=40)
  par(mar=c(2.1,1.6,4.1,0.1), mfrow=c(1,2), oma=c(7,0,2,0), xpd=NA, bty="n")
  plot(pot_eg, col=pn_cols, xaxt="n", yaxt="n", legend=FALSE)
  axis(1, at=c(147.45, 147.55), labels=c("",""), lwd=3)
  mtext("0.1 degrees", 1, 1, cex=1.3)
  plot(maskras, col=alpha("orange", 0.3), add=TRUE, legend=FALSE)
  # points(mid[1], mid[2], pch=4, lwd=5, cex=1.5)
  points(mid[1], mid[2], cex=1.3, col="orange", pch=15)
  
  plot(guelphia_potential, col=pn_cols, xaxt="n", yaxt="n", legend=FALSE)
  plot(coarse_catch, col=alpha("orange", 0.3), add=TRUE, legend=FALSE)
  points(allxy[45, 1], allxy[45,2], col="orange", pch=15, cex=4)
  axis(1, at=c(-37.12, -37.02), labels=c("",""), lwd=3)
  mtext("0.1 degrees", 1, 1, cex=1.3)
  
  mtext("Catchment around a point", 3, outer=TRUE, line=-1, font=2, cex=1.6)
  par(mfrow=c(1,1), new=TRUE, mar=c(0,0,0,0), oma=c(0,0,0,0))
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  legend(0.8,0.15, fill=c("orange", alpha("orange", 0.3)), c("Site","Pixels in \ncatchment"), 
         cex=1.4, bty="n")
  plot(guelphia_potential, col=pn_cols, legend.only=TRUE,
       axis.args=list(at=c(minValue(guelphia_potential), maxValue(guelphia_potential)), 
                      labels=c("Low","High"), cex.axis=1.4),
       legend.args=list("JEV transmission suitability", side=3, line=1, cex=1.4),
       legend.mar=3,
       legend.shrink=0.3,
       horizontal=TRUE)
  dev.off()}



plot(guelphia_potential, col=pn_cols, legend.only=TRUE,
     axis.args=list(at=c(minValue(pot_eg), maxValue(pot_eg)), 
                    labels=c("Low","High"), cex.axis=1.4),
     legend.args=list("JEV transmission suitability", side=2, line=2, cex=1.4))









