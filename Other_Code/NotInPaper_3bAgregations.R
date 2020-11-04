##################
#predmat <- matrix(5/42, 7, 7)
#diag(predmat) <- 2/7
#colnames(predmat) <- spec_names
#newdata <- list(mean_comp = rbind(1/7, predmat), indicator = rep(1,8))
#
#res <- predict(metamod, newdata)
#RR <- res[-1] / res[1]
#
#plot(RR, -(1:7), pch = 16, col = spec_pal)
#abline(v = 1)
#
#plot(exp(coef_all), -(1:7), pch = 16, col = spec_pal)
#abline(v = 1)

#-------------------------------------
#  Aggregated composition NOT IN PAPER
#-------------------------------------
# The aggregation with consider
agg_inds <- list(INOR = 1:3, CARBON = 4:5, SS = 6, DUST = 7)
pa <- length(agg_inds)
agg_pal <- c("purple", "darkgreen", 5, 7)

# Create the aggregation
agg_comp <- t(sapply(dlist_spec, function(x){
    x <- x[,spec_inds]
    # Aggregating
    xs <- sapply(agg_inds, function(inds) rowSums(x[,inds, drop = F]))
    # Zero value imputation as in Martín-Fernández et al. (2003)
    if (any(xs == 0)) xs <- multRepl(xs, label = 0, dl = rep(1e-5, pa))
    # Compositional mean (pass through the irl transformation)
    mean(acomp(xs))
}))

# We consider the component "NATURAL" to be the baseline and thus estimated the effect of two other compared to this one
#   This leads to consider the alr of Aitchison 
#   Note that this reduces to considering the approach of Aitchison & Bacon-Shone (1984) 

agg_reg <- mixmeta(coefall ~ alr(agg_comp) + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)

coef_agg <- coef(agg_reg)[2:4]
se_agg <- sqrt(diag(vcov(agg_reg))[2:4])
lo_agg <- coef_agg - 1.96 * se_agg
up_agg <- coef_agg + 1.96 * se_agg
sig_agg <- lo_agg > 0 | up_agg < 0

# Forestplot

x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(coef_agg, -seq_len(pa - 1), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(c(lo_agg, up_agg)))
abline(v = 0)
segments(lo_agg, -seq_len(pa - 1), up_agg, -seq_len(pa - 1), 
  col = "darkgrey", lwd = 2)
points(coef_agg, -seq_len(pa - 1), cex = ifelse(sig_agg, 1.5, 1), 
   pch = ifelse(sig_agg, 15, 16), col = agg_pal[-pa])
axis(1)
axis(2, at = -seq_len(pa - 1), labels = names(agg_inds)[-pa], 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)

dev.print(png, filename = "Results/3f_forestplot_aggregated.png", 
  units = "in", res = 100)