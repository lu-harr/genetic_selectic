gradient_legend = function(cuts, cols, xticks=c(), leg_pos="right"){
  # doesn't quite work quite right yet
  lgd_ = rep(NA, length(levels(cuts)))
  if (length(xticks) == 0){
    pick_out = seq(1, length(levels(cuts)), length.out = 5)
    to_label = levels(cuts)[pick_out]
    to_label = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", to_label)) %>%
      signif(digits = 3)# manipulate char from cut()
    lgd_[pick_out] = to_label
  }
  legend(leg_pos,
         legend = rev(lgd_),
         fill = rev(cols),
         border = NA,
         y.intersp = 0.5)
}

pr_chop = function(x){
  return(ifelse(x > 1, 1, x))
}

# how does objective change with
ras = raster(nrow=500, ncol=500, xmn=0, xmx=1, ymn=0, ymx=1)
values(ras) = 0

# fiddling with products/sums?
# rass = stack(ras, ras, ras, ras, ras)
# names(rass) = c("curr","prop","chop_prop", "dont_chop_prop", "two_chopped_props")
# locs = xyFromCell(ras, 1:ncell(ras))
# values(rass$curr) = locs[,1] # changes across xax
# values(rass$prop) = locs[,2] # changes across yax
# values(rass$chop_prop) = apply(cbind(values(rass$curr),
#                                      values(rass$prop)),
#                                1, function(x){obj_to_pr(x)})
# values(rass$dont_chop_prop) = values(rass$prop) / values(rass$curr)
# 
# values(rass$two_chopped_props) = values(rass$chop_prop)**2
# 
# par(mfrow=c(1,2))
# plot(rass$two_chopped_props)
# plot(rass$chop_prop)

# still fiddling with products/sums ...
# obj1 = seq(0,1,length.out=500)
# plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="Obj1", ylab="Utility = obj1*obj2*obj3")
# sapply(seq(0,1,length.out=10),
#        function(x){lines(obj1, obj1*x)}) # obj3 = 1
# sapply(seq(0,1,length.out=10),
#        function(x){lines(obj1, obj1*x*0.9, col="blue")}) # obj3 = 0.9
# sapply(seq(0,1,length.out=10),
#        function(x){lines(obj1, obj1*x*0.8, col="red")}) # obj3 = 0.8
# When total utility is 0.5, what are the values of the objs?
# (0.5, 1), (0.5555556, 0.9)

# create some example objs
noisy = data.frame(obj1 = rexp(1000),
                   obj2 = rexp(1000),
                   obj3 = rexp(1000))

{png("~/Desktop/knowlesi/multi_site/exploratory_plots/truncate_intersect.png",
     height=1200, 
     width=2000, 
     pointsize=25)
  par(mfrow=c(1,2), oma=c(0,0,2,0))

# highlight points good for both (INTERSECT)
best_points = which(noisy$obj1 > 0.8 & noisy$obj2 > 0.8)
okay_points = which(noisy$obj1 > 0.5 & noisy$obj2 > 0.5)

# apply probability chop before product
plot(noisy$obj1,
     pr_chop(noisy$obj1) * pr_chop(noisy$obj2),
     xlab="First objective", ylab="Utility = chopped(first obj) * chopped(second obj)",
     main="Apply truncation before\n product operation")
points(noisy$obj1[okay_points],
     pr_chop(noisy$obj1[okay_points]) * pr_chop(noisy$obj2[okay_points]), col="blue", lwd=2)
points(noisy$obj1[best_points],
       pr_chop(noisy$obj1[best_points]) * pr_chop(noisy$obj2[best_points]), col="red", lwd=2)
abline(a=0, b=1)

# apply probability chop after product
plot(noisy$obj1,
     pr_chop(noisy$obj1 * noisy$obj2),
     xlab="First objective", ylab="Utility = chopped(first obj * second obj)",
     main="Apply truncation after\n product operation")
points(noisy$obj1[okay_points],
       pr_chop(noisy$obj1[okay_points] * noisy$obj2[okay_points]), col="blue", lwd=2)
points(noisy$obj1[best_points],
       pr_chop(noisy$obj1[best_points] * noisy$obj2[best_points]), col="red", lwd=2)
abline(a=0, b=1)
legend("bottomright",
       fill = c("red", "blue"),
       c("obj1 > 0.8 & obj2 > 0.8", "obj1 > 0.5 & obj2 > 0.5"))
mtext("How does truncation step affect alternatives good for *both* objectives?", outer = TRUE)

dev.off()}
# so the chopping process definitely demotes alternatives that are good for *both* objectives ... how about either?

# need to think some more about this one
plot(pr_chop(noisy$obj1) * pr_chop(noisy$obj2),
     pr_chop(noisy$obj1 * noisy$obj2),
     xlab="Utility = chopped(first obj) * chopped(second obj)",
     ylab="Utility = chopped(first obj * second obj)",
     main="Utility same or greater when \ntruncation applied after product",
     col=viridis(100)[cut(log(noisy$obj1), breaks=100)], pch=16)


{png("~/Desktop/knowlesi/multi_site/exploratory_plots/truncate_union.png",
     height=1200, 
     width=2000, 
     pointsize=25)
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  
best_points = which(noisy$obj1 > 0.8 | noisy$obj2 > 0.8)
okay_points = which(noisy$obj1 > 0.5 | noisy$obj2 > 0.5)

# a bit more continuous
min_obj = sapply(1:nrow(noisy),
                 function(x){
                   min(noisy$obj1[x], noisy$obj2[x])
                 })

# apply probability chop before product
plot(noisy$obj1,
     pr_chop(noisy$obj1) * pr_chop(noisy$obj2),
     xlab="Factor 1", ylab="Alpha = truncated(factor 1) * chopped(factor 2)",
     main="Apply truncation before\n product operation",
     col=viridis(100)[cut_number(min_obj, n=100)],
     pch=16)
# points(noisy$obj1[okay_points],
#        pr_chop(noisy$obj1[okay_points]) * pr_chop(noisy$obj2[okay_points]), col="blue")
# points(noisy$obj1[best_points],
#        pr_chop(noisy$obj1[best_points]) * pr_chop(noisy$obj2[best_points]), col="red")
points(noisy$obj1,
       pr_chop(noisy$obj1) * pr_chop(noisy$obj2), 
       col=alpha("red", noisy$obj2/max(noisy$obj2)), lwd=4, cex=1.5)
abline(a=0, b=1)

# apply probability chop after product
plot(noisy$obj1,
     pr_chop(noisy$obj1 * noisy$obj2),
     xlab="Factor 1", ylab="Alpha = truncated(factor 1 * factor 2)",
     main="Apply truncation after\n product operation",
     col=viridis(100)[cut_number(min_obj, n=100)],
     pch=16)
# points(noisy$obj1[okay_points],
#        pr_chop(noisy$obj1[okay_points] * noisy$obj2[okay_points]), col="blue")
# points(noisy$obj1[best_points],
#        pr_chop(noisy$obj1[best_points] * noisy$obj2[best_points]), col="red")
points(noisy$obj1,
       pr_chop(noisy$obj1 * noisy$obj2), 
       col=alpha("red", noisy$obj2/max(noisy$obj2)), lwd=4, cex=1.5)
# points(noisy$obj1,
#        pr_chop(noisy$obj1 * noisy$obj2), 
#        col=alpha("blue", noisy$obj1/max(noisy$obj1)), lwd=4, cex=2)
gradient_legend(cut_number(min_obj, n=20), cols = viridis(20))
mtext("min(factor 1, factor 2)", side=4)
abline(a=0, b=1)
mtext("How does chopping affect alternatives good for *one* of the objectives?", outer = TRUE)

dev.off()}

##########################################################################
# Here's another one ... raster with contour

ras = raster(nrow=500, ncol=500, xmn=0, xmx=2, ymn=0, ymx=2)
locs = xyFromCell(ras, 1:ncell(ras))
obj1 = locs[,1]
obj2 = locs[,2]

{png("~/Desktop/knowlesi/multi_site/exploratory_plots/truncate_raster.png",
     height=1000, 
     width=2000, 
     pointsize=25)
  par(mfrow=c(1,2))

par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 5.1))
values(ras) = pr_chop(obj1 * obj2)
#values(ras) = obj1*obj2
plot(ras, main="Acceptance probability: Post-truncate", xlab="Factor 1", 
     ylab="Factor 2", col=viridis(100))
abline(h = 1, lty = 2)
abline(v = 1, lty = 2)
lines(c(0.5, 0.5, seq(0.5,1,length.out=50), 1, 2),
      c(2, 1, 0.5/seq(0.5,1,length.out=50), 0.5, 0.5),
      col="grey60")
lines(seq(0.1, 2, length.out=100),
      0.5/seq(0.1, 2, length.out=100), col="white", lwd=2)


values(ras) = pr_chop(obj1) * pr_chop(obj2)
plot(ras, main="Acceptance probability: Pre-truncate", xlab="Factor 1", 
     ylab="Factor 2", col=viridis(100))
abline(h = 1, lty = 2)
abline(v = 1, lty = 2)
lines(seq(0.1, 2, length.out=100),
      0.5/seq(0.1, 2, length.out=100),
      col="grey60")
lines(c(0.5, 0.5, seq(0.5,1,length.out=50), 1, 2),
      c(2, 1, 0.5/seq(0.5,1,length.out=50), 0.5, 0.5), col="white", lwd=2)
dev.off()}
# interesting ...


##########################################################################













