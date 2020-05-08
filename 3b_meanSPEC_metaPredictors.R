#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#                 Meta-Predictors PM2.5 SPEC 
#
#####################################################################

library(mixmeta)

setwd("C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Projects/MCC-HetPoll")
load("Data/2_FirstStageResults.RData")

#-------------------------------------
#   Preparing constituent data
#-------------------------------------
spec_inds <- grep("PM25", colnames(dlist_spec[[1]]))
spec_names <- sapply(strsplit(colnames(dlist_spec[[1]])[spec_inds], "_"),
  "[", 2)
spec_pal <- c("red", "pink", "blue", "black", "green", "cyan", "yellow")

mean_spec <- t(sapply(dlist_spec, function(d){
  colMeans(d[,spec_inds])
}))

#-------------------------------------
#   Meta-regression model
#-------------------------------------

metamod <- mixmeta(coefall ~ mean_spec, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml")

beta_spec <- coef(metamod)[-1]
se_spec <- sqrt(diag(vcov(metamod))[-1])
lo_spec <- beta_spec - 1.96 * se_spec
up_spec <- beta_spec + 1.96 * se_spec
sig_spec <- lo_spec > 0 | up_spec < 0
bord <- order(beta_spec, decreasing = T)

# Forestplot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(beta_spec[bord], -seq_along(beta_spec), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(c(lo_spec, up_spec)))
segments(lo_spec[bord], -seq_along(beta_spec), up_spec[bord],
  -seq_along(beta_spec), col = "darkgrey", lwd = 2)
points(beta_spec[bord], -seq_along(beta_spec), 
   col = spec_pal[bord], cex = ifelse(sig_spec[bord], 1.5, 1), 
   pch = ifelse(sig_spec[bord], 15, 16))
axis(1)
axis(2, at = -seq_along(beta_spec), labels = spec_names[bord], 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3b_forestplot_betaSPEC.png", 
  units = "in", res = 100)