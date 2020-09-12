#####################################################################
#
#                 MCC HetPoll
#             Bonus: supplementary plots and results
#
#####################################################################

library(Ternary)
library(compositions)
library(mixmeta)
library(fields)
library(colorspace)
library(fda)

load("Data/3_MetaResults.RData")

#-------------------------------------
#  Ternary plot
#-------------------------------------

# Parameters
tern_comp <- c(1, 3, 4) # Which components
res <- 40 # resolution of the plot

colpal <- colorRampPalette(c("blue", "white", "red"))

# Create function to predict for Ternary plot
tern_pred <- function(c1, c2, c3, tern_comp, p, ov_mean, metamod, cols){
  # Reduce the chosen components to their average composition
  cgrid <- cbind(c1, c2, c3) * sum(ov_mean[tern_comp])
  
  # Initiate data  
  xc <- matrix(NA, nrow = nrow(cgrid), p, dimnames = list(NULL, spec_names))
  xc[,tern_comp] <- data.matrix(cgrid)
    
  # Add other components overall mean
  xc[,-tern_comp] <- t(replicate(nrow(cgrid), ov_mean[-tern_comp]))
    
  # Predictions
  newdat <- data.frame(indicator = rep(0, nrow(cgrid)))
  newdat$alr_comp <- alr(xc)
  tern_preds <- exp(predict(metamod, newdat))
  
  # Create color  
  colpred <- cols[cut(tern_preds, 
    breaks = seq(0.985, 1.015, length.out = length(cols) + 1))
  ]
  attr(colpred, "predictions") <- tern_preds
  colpred
}

tern_values <- TernaryPointValues(tern_pred, direction = 1,
  tern_comp = tern_comp, p = p, ov_mean = ov_mean, metamod = metamod,
  cols = colpal(20), resolution = 100)

# Function for observed compositions
crngs <- apply(mean_comp[,tern_comp], 2, range)
obs_func <- function(a, b, c){
  is_obs <- a >= crngs[1, 1] & a <= crngs[2, 1] &
    b >= crngs[1, 2] & b <= crngs[2, 2] &
    c >= crngs[1, 3] & c <= crngs[2, 3]
  as.numeric(is_obs)
}


#----- Plot ternary diagram with predictions -----
x11()
par(mar = c(5, 4, 4, 7) + .1)
TernaryPlot(alab = spec_labs[tern_comp[1]], blab = spec_labs[tern_comp[2]],
  clab = spec_labs[tern_comp[3]], lab.col = spec_pal[tern_comp],
  grid.col = "white", grid.minor.col = "white")
ColourTernary(tern_values, spectrum = NULL)
TernaryPoints(mean_comp[,tern_comp], pch = 16)

# Adding the color scale
image.plot(zlim = c(0.985, 1.015), col = colpal(20), 
  legend.only = T)
mtext("RR", side = 3, at = par("usr")[2] + .15, xpd = T, cex = 1.5, line = -1)

dev.print(png, filename = sprintf("Paper_Figures/SuppFigure1_TernaryPredictions_%s.png", 
    paste(spec_names[tern_comp], collapse = "_")), 
  units = "in", res = 100)
dev.print(pdf, file = sprintf("Paper_Figures/SuppFigure1_TernaryPredictions_%s.pdf", 
    paste(spec_names[tern_comp], collapse = "_")))

#-------------------------------------
#  Forest plot of the estimated coefficients in the meta-regression model
#-------------------------------------

# Coefficient retrieving
coef_est <- rep(NA, p)
coef_est[-p] <- coef(metamod)[2:p]
coef_est[p] <- -sum(coef_est, na.rm = T) # beta7 = - sum(beta1, ..., beta6)

# Standard errors
se_est <- rep(NA, p)
se_est[-p] <- sqrt(diag(vcov(metamod))[2:p])
se_est[p] <- sqrt(sum(vcov(metamod)[2:p,2:p])) # var(beta7) = sum_j sum_k cov(betaj, betak) (variance of sum of random variables)

# Confidence intervals
lo_est <- coef_est - 1.96 * se_est
up_est <- coef_est + 1.96 * se_est
sig_est <- lo_est > 0 | up_est < 0

# We mulitply everything by ln(2) to interpret it as the impact of doubling the relative proportion of the composnents
coef_est <- log(2) * coef_est
lo_est <- log(2) * lo_est
up_est <- log(2) * up_est

# Forestplot of the RR increase
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(exp(coef_est), -seq_len(7), 
  xlab = "Relative Excess Risk", ylab = "", cex.lab = 1.3,
  axes = F, xlim = range(c(exp(lo_est), exp(up_est))))
abline(v = 1)
segments(exp(lo_est), -seq_len(7), exp(up_est), -seq_len(7), lwd = 2,
  col = grey(.2))
points(exp(coef_est), -seq_len(7), cex = ifelse(sig_est, 2, 1.5), 
   pch = ifelse(sig_est, 15, 16), col = spec_pal)
axis(1, cex.axis = 1.1)
axis(2, at = -seq_len(7), labels = spec_labs, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0, cex.axis = 1.3)

dev.print(png, filename = "Paper_Figures/SuppFigure2_ForestPlot.png", 
  units = "in", res = 200)
dev.print(pdf, file = "Paper_Figures/SuppFigure2_ForestPlot.pdf")


#-------------------------------------
# Checking meta-model assumptions
#-------------------------------------

#---- Normality of residuals ----
x11()
hist(residuals(metamod), border = "white", col = 4, xlab = "Residuals", freq = F,
  main = "")

dev.print(png, filename = "Paper_Figures/SuppFigure3_ResidHist.png", 
  units = "in", res = 200)
dev.print(pdf, file = "Paper_Figures/SuppFigure3_ResidHist.pdf")

#---- Plot residuals by region ----

# Organise by regions
nregion <- length(unique(cities$region))
reg_ord <- order(cities$region)
cities_ord <- cities[reg_ord,]
cities_ord$region <- droplevels(cities_ord$region)
res_ord <- residuals(metamod)[reg_ord]

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
  
dev.print(png, filename = "Paper_Figures/SuppFigure4_ResidRegion.png", 
  units = "in", res = 200)
dev.print(pdf, file = "Paper_Figures/SuppFigure4_ResidRegion.pdf")

#---- Plot residuals by components ----
x11(height = 10, width = 10)
layout(pm)
for (j in seq(p)){
  plot(100 * mean_comp[,j], residuals(metamod), col = spec_pal[j], pch = 16, 
    xlim = 100 * range(mean_comp), ylim = range(residuals(metamod)),
    ylab = "Residuals", xlab = "Proportion (%)",
    main = bquote(bold(.(spec_labs[j][[1]]))))
}

dev.print(png, filename = "Paper_Figures/SuppFigure5_ResidComponent.png", 
  units = "in", res = 200)
dev.print(pdf, file = "Paper_Figures/SuppFigure5_ResidComponent.pdf")

#---- Check extreme residuals ----
res_sort <- order(residuals(metamod))
# Negative residuals
cities[res_sort[1:3],]
# Positive residuals
cities[res_sort[nrow(cities) - 0:2],]