#####################################################################
#
#                 MCC HetPoll
#             Bonus: supplementary plots and results
#
#####################################################################

library(Ternary)
library(compositions)
library(zCompositions)
library(mixmeta)
library(fields)
library(colorspace)
library(fda)

load("Data/3_MetaResults.RData")

#-------------------------------------
#  Supplemental Material A: Indicators PCA
#-------------------------------------

var_prop <- 100 * summary(pca_indic)$importance[2,]

#----- Screeplot -----
x11(width = 15)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.2, mar = c(5, 4, 5, 2) + .1)
plot(var_prop, type = "b", 
  xlab = "Principal component", ylab = "Variance proportion (%)",
  pch = 16, xaxt = "n", main = "Scree plot")
axis(1, at = 1:length(var_prop))

#----- Biplot -----
# Color by MCC region
nregion <- length(unique(cities$region))
reg_pal <- rainbow_hcl(nregion)
reg_pch <- rep_len(15:18, nregion)
names(reg_pal) <- names(reg_pch) <- unique(cities$region)

# Loadings
loads <- pca_indic$rotation[,1:2] * sum(pca_indic$sdev)

plot(pca_indic$x[,1], pca_indic$x[,2], asp = 1, 
  xlab = sprintf("PC 1 (%2.1f %%)", var_prop[1]), 
  ylab = sprintf("PC 2 (%2.1f %%)", var_prop[2]),
  pch = reg_pch[as.character(cities$region)], 
  col = reg_pal[as.character(cities$region)],
  main = "Biplot"
)
abline(v = 0, lty = 2, col = "darkgrey")
abline(h = 0, lty = 2, col = "darkgrey")
arrows(0, 0, loads[,1], loads[,2], length = .1, lwd = 2)
text(loads[,1], loads[,2], rownames(loads), pos = 2 + sign(loads[,2]))
legend("topright", legend = unique(cities$region), 
  pch = reg_pch[as.character(unique(cities$region))], 
  col = reg_pal[as.character(unique(cities$region))],
  inset = .02, ncol = 1, xpd = T
)

dev.print(png, filename = "Results/FigureS1.png", 
  units = "in", res = 100)
dev.print(pdf, file = "Results/FigureS1.pdf")

#-------------------------------------
#  Supplemental Material B: Residual Analysis
#-------------------------------------

resvec <- rep(NA, nrow(cities))
resvec[conv] <- residuals(metamod)

#---- Normality of residuals ----
x11()
hist(resvec, border = "white", col = 4, xlab = "Residuals", freq = F,
     main = "")

dev.print(png, filename = "Results/FigureS2.png", 
          units = "in", res = 200)
dev.print(pdf, file = "Results/FigureS2.pdf")

#---- Plot residuals by region ----

# Organise by regions
nregion <- length(unique(cities$region))
reg_ord <- order(cities$region)
cities_ord <- cities[reg_ord,]
cities_ord$region <- droplevels(cities_ord$region)
res_ord <- resvec[reg_ord]

# Color palette and point type
reg_pal <- rainbow_hcl(nregion)
reg_pch <- rep_len(15:18, nregion)

# limits of each region
change_reg <- which(diff(as.numeric(cities_ord$region)) != 0) + .5

# plot residuals
x11(width = 10)
par(mar = c(5, 4, 4, 10))
plot(res_ord, pch = reg_pch[cities_ord$region], col = reg_pal[cities_ord$region],
     xlab = "", ylab = "Residual")
abline(h = 0)
legend(par("usr")[2], par("usr")[4], unique(cities_ord$region), pch = reg_pch,
       col = reg_pal, bty = "n", xpd = T, title = "Region")

dev.print(png, filename = "Results/FigureS3.png", 
          units = "in", res = 200)
dev.print(pdf, file = "Results/FigureS3.pdf")

#---- Plot residuals by components ----
pm <- matrix(c(1:6, 0, 7, 0), nrow = 3, byrow = T)

x11(height = 10, width = 10)
layout(pm)
for (j in seq(p)){
  plot(100 * mean_comp[,j], resvec, col = spec_pal[j], pch = 16, 
       xlim = 100 * range(mean_comp), ylim = range(resvec, na.rm = T),
       ylab = "Residuals", xlab = "Proportion (%)",
       main = bquote(bold(.(spec_labs[j][[1]]))))
}

dev.print(png, filename = "Results/FigureS4.png", 
          units = "in", res = 200)
dev.print(pdf, file = "Results/FigureS4.pdf")

#---- Check extreme residuals ----
res_sort <- order(resvec)
# Negative residuals
cities[res_sort[1:3],]
# Positive residuals
cities[res_sort[nrow(cities) - 0:2],]


#-------------------------------------
#  Supplemental Material C: Alternative representation of results
#-------------------------------------

#----- Figure S5: Predictions on the same panel -----

x11(width = 10, height = 7)
par(mar = c(5, 4, 4, 10) + .1)
matplot(100 * cseq, preds, type = "l", lty = 2, col = spec_pal, 
  ylab = "RR", xlab = "Proportion (%)", cex.lab = 1.2)
# Add thicker lines for the observed range
matlines(100 * cseq, preds_obs, lty = 1, lwd = 3, col = spec_pal)
# Add legends
abline(h = 1)
lg <- legend(par("usr")[2], par("usr")[4], spec_labs, lwd = 3, col = spec_pal, 
  bty = "n", xpd = T, y.intersp = 1.2)
legend(par("usr")[2], with(lg$rect, top - h), c("Observed range", "Extrapolation"), 
  lwd = c(3, 1), lty = 1:2, bty = "n", xpd = T, y.intersp = 1.2, col = grey(.5))
  
dev.print(png, filename = "Results/FigureS5.png", 
  units = "in", res = 100)
dev.print(pdf, file = "Results/FigureS5.pdf")

#----- Plot ternary diagram with predictions -----
## Create function to predict for Ternary plot
tern_pred <- function(c1, c2, c3, tern_comp, p, ov_mean, metamod, cols){

  # Reduce the chosen components to their average composition
  cgrid <- cbind(c1, c2, c3) * sum(ov_mean[tern_comp])
  
  # Initiate data  
  xc <- matrix(NA, nrow = nrow(cgrid), p, dimnames = list(NULL, spec_names))
  xc[,tern_comp] <- data.matrix(cgrid)
    
  # Add other components overall mean
  xc[,-tern_comp] <- t(replicate(nrow(cgrid), ov_mean[-tern_comp]))
    
  # Predictions
  newdat <- list(indicators = matrix(0, nrow(xc), 2, 
    dimnames = list(NULL, sprintf("PC%i", 1:2))))
  newdat$alr_comp <- alr(xc)
  tern_preds <- exp(predict(metamod, newdat))
  
  # Create color  
  colpred <- cols[cut(tern_preds, 
    breaks = seq(0.985, 1.015, length.out = length(cols) + 1))
  ]
  attr(colpred, "predictions") <- tern_preds
  colpred
}

## Parameters
res <- 40 # resolution of the plot
colpal <- colorRampPalette(c("blue", "white", "red"))

tern_comps <- list(c(2,3,7), 1:3, 4:6, 5:7, c(1,3,4))  # Which components

## create the plot
x11(height = 15, width = 10)
par(mfrow = n2mfrow(length(tern_comps) + 1, asp = .6), mar = c(1, 1, 1, 1) + .1)
for (i in seq_along(tern_comps)){
  tern_values <- TernaryPointValues(tern_pred, direction = 1,
    tern_comp = tern_comps[[i]], p = p, ov_mean = ov_mean, metamod = metamod,
    cols = colpal(20), resolution = 100)
  
  # Function for observed compositions
  crngs <- apply(mean_comp[,tern_comps[[i]]], 2, range)
  obs_func <- function(a, b, c){
    is_obs <- a >= crngs[1, 1] & a <= crngs[2, 1] &
      b >= crngs[1, 2] & b <= crngs[2, 2] &
      c >= crngs[1, 3] & c <= crngs[2, 3]
    as.numeric(is_obs)
  }
  
  # Add to plot
  TernaryPlot(alab = spec_labs[tern_comps[[i]][1]], blab = spec_labs[tern_comps[[i]][2]],
    clab = spec_labs[tern_comps[[i]][3]], lab.col = spec_pal[tern_comps[[i]]],
    grid.col = "white", grid.minor.col = "white")
  ColourTernary(tern_values, spectrum = NULL)
  TernaryPoints(mean_comp[,tern_comps[[i]]], pch = 16)
}

# Adding the color scale
plot.new()
image.plot(zlim = c(0.985, 1.015), col = colpal(20), 
  legend.only = T, xpd = T, smallplot = c(0, .05, 0.2, .8))
text(0, .92, "RR", xpd = T, cex = 1.5)

dev.print(png, filename = "Results/FigureS6.png", units = "in", res = 100)
dev.print(pdf, file = "Results/FigureS6.pdf")


#-------------------------------------
# Table S2 : Correlation table between total pm and absolute components
#-------------------------------------

# Correlation between absolute components and total_pm
tot_spec <- do.call(rbind, dlist_spec)[,spec_inds]
cormat <- cor(rowSums(tot_spec), tot_spec)

# Correlation between relative components and mean_pm
cormat <- rbind(cormat, cor(mean_pm, mean_comp))

# Export
rownames(cormat) <- c("Absolute", "Relative")
colnames(cormat) <- spec_names
cormat <- round(cormat, 2)

write.table(cormat, file = "Results/TableS2.csv", quote = F,
  row.names = T, sep = ";")