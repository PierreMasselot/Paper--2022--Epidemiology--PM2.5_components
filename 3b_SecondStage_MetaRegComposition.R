#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#              experiment with mixtures
#
#####################################################################

# Here I simply consider the approach of Aitchison & Bacon-Shone (1984)
# It consists in transforming the constituents with ALR and using them in a regression model.


library(mixmeta)
library(compositions)
library(zCompositions)

load("Data/2_FirstStageResults.RData")

#-------------------------------------
#   Preparing constituent objects
#-------------------------------------

spec_inds <- grep("PM25", colnames(dlist_spec[[1]]))
spec_names <- c("SO4", "NH4", "NO3", "BC", "OC", "SS", "DUST")
spec_labs <- c(expression(SO[4]^{"2-"}), expression(NH[4]^{"+"}), expression(NO[3]^{"-"}), 
  "BC", "OC", "SS", "DUST")
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# Mean per city
mean_comp <- sapply(dlist_spec, function(x){
  # Zero value imputation as in Martín-Fernández et al. (2003)
  imp <- multRepl(x[,spec_inds], label = 0, dl = rep(1e-5, 7))
  # Transformation to compositional object
  xc <- acomp(imp)
  # Compositional mean (pass through the irl transformation)
  mean(xc)
})
mean_comp <- t(mean_comp)
p <- ncol(mean_comp)
colnames(mean_comp) <- spec_names

# Transforming as ALR
alr_comp <- alr(mean_comp)

#-------------------------------------
#   Scores of MCC indicators
#-------------------------------------

### Which indicators
# Density has too much NAs
# For the UCD indicators, we take both the 2000 and 2014/2015 ones
indic_names <- c("oldpopprop", "GDP", "Poverty", "avgtmean", 
  "totalrange", "E_GR_AV00", "E_GR_AV14", "B00", "B15")
indic_form <- sprintf("~ %s", paste(indic_names, collapse = " + "))

# PCA
pca_indic <- princomp(as.formula(indic_form), data = mcc.indicators,
  na.action = na.exclude)
# screeplot(pca_indic)
# The first component accounts for 99% of the variance
indicator <- pca_indic$scores[,1]

#-------------------------------------
#  All components
#-------------------------------------

# Model with ALR as meta-predictors
metamod <- mixmeta(coefall ~ alr_comp + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
  
summary(metamod) # To obtain Q and I2 statistics

# Coefficient retrieving
coef_all <- rep(NA, p)
coef_all[-p] <- coef(metamod)[2:p]
coef_all[p] <- -sum(coef_all, na.rm = T) # beta7 = - sum(beta1, ..., beta6)

# Standard errors
se_all <- rep(NA, p)
se_all[-p] <- sqrt(diag(vcov(metamod))[2:p])
se_all[p] <- sqrt(sum(vcov(metamod)[2:p,2:p])) # var(beta7) = sum_j sum_k cov(betaj, betak) (variance of sum of random variables)

# Confidence intervals
lo_all <- coef_all - 1.96 * se_all
up_all <- coef_all + 1.96 * se_all
sig_all <- lo_all > 0 | up_all < 0

# We mulitply everything by ln(2) to interpret it as the impact of doubling the relative proportion of the composnents
coef_all <- log(2) * coef_all
lo_all <- log(2) * lo_all
up_all <- log(2) * up_all

# Forestplot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(coef_all, -seq_len(7), 
  xlab = "Log(RR) increase", ylab = "", cex.lab = 1.3,
  axes = F, xlim = range(c(lo_all, up_all)))
abline(v = 0)
segments(lo_all, -seq_len(7), up_all, -seq_len(7), lwd = 2,
  col = grey(.2))
points(coef_all, -seq_len(7), cex = ifelse(sig_all, 2, 1.5), 
   pch = ifelse(sig_all, 15, 16), col = spec_pal)
axis(1, cex.axis = 1.1)
axis(2, at = -seq_len(7), labels = spec_labs, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0, cex.axis = 1.3)

dev.print(png, filename = "Paper_Figures/Figure3.png", 
  units = "in", res = 200)
  
# Alternative forestplot: RR increase
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(exp(coef_all), -seq_len(7), 
  xlab = "Relative Excess Risk", ylab = "", cex.lab = 1.3,
  axes = F, xlim = range(c(exp(lo_all), exp(up_all))))
abline(v = 1)
segments(exp(lo_all), -seq_len(7), exp(up_all), -seq_len(7), lwd = 2,
  col = grey(.2))
points(exp(coef_all), -seq_len(7), cex = ifelse(sig_all, 2, 1.5), 
   pch = ifelse(sig_all, 15, 16), col = spec_pal)
axis(1, cex.axis = 1.1)
axis(2, at = -seq_len(7), labels = spec_labs, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0, cex.axis = 1.3)

dev.print(png, filename = "Paper_Figures/Figure3_RR.png", 
  units = "in", res = 200)
  

#-------------------------------------
#  Prediction of the RR for single components
#-------------------------------------

# We create a sequence for the component
cseq <- seq(.01, .99, by = .01)

# Overall mean of composition
ov_mean <- mean(acomp(mean_comp))

# Prepare prediction data.frame
newdat <- data.frame(indicator = rep(0, length(cseq)))

# Prepare objects to store predictions and confidence limits
preds <- plo <- pup <- matrix(NA, length(cseq), p)
for (j in seq_len(p)){
  # Create a compositional grid. The component of interests varies from 0 to 1 (excluded) and the other are taken as the overall mean adjusted for closure while keeping the subcomposition constant
  xj <- matrix(NA, length(cseq), p, dimnames = list(NULL, spec_names))
  xj[,-j] <- sapply(ov_mean[-j] / sum(ov_mean[-j]), "*", 1 - cseq)
  xj[,j] <- cseq
  
  # Transform by ALR and add to the newdat data.frame
  newdat$alr_comp <- alr(xj)
  
  # Make predictions
  pred_obj <- predict(metamod, newdata = newdat, ci = T)
  
  # Store predictions and prediction intervals only in the range of observed values
  obs_rng <- range(mean_comp[,j])
  is_obs <- cseq >= obs_rng[1] & cseq <= obs_rng[2]
  
  preds[is_obs,j] <- exp(pred_obj[is_obs,1])
  plo[is_obs,j] <- exp(pred_obj[is_obs,2])
  pup[is_obs,j] <- exp(pred_obj[is_obs,3])
}

#--- Plot all predictions in different panels ----
# Panel matrix
pm <- matrix(rep(c(1:6, 0, 7, 0), each = 2), nrow = 3, byrow = T)

# Plot
x11(height = 10, width = 10)
par(mar = c(5, 4, 3, 1), cex.main = 1.5, cex.lab = 1.2)
layout(pm)
for (j in seq_len(p)){
  # Initiate an empty plot with labels
  plot(0, 0, col = NA, xlim = 100 * range(mean_comp), 
    ylim = c(min(plo, na.rm = T), max(pup, na.rm = T)),
    ylab = "RR", xlab = "Proportion (%)",
    main = bquote(bold(.(spec_labs[j][[1]]))))
  # Draw prediction interval
  nona <- !is.na(preds[,j])
  polygon(100 * c(cseq[nona], rev(cseq[nona])), c(plo[nona,j], rev(pup[nona,j])), 
    col = adjustcolor(spec_pal[j], .2), border = NA)
  # Add prediction
  lines(100 * cseq, preds[,j], lty = 1, lwd = 2, col = spec_pal[j])
  abline(h = 1)
}

dev.print(png, filename = "Results/3_RRpredictions.png", 
  units = "in", res = 100)
dev.print(pdf, file = "Results/3_RRpredictions.pdf")

#-------------------------------------
#  Ternary plot
#-------------------------------------

# Parameters
tern_comp <- c(2, 3, 7) # Which components
res <- 40 # resolution of the plot

colpal <- colorRampPalette(c("blue", "white", "red"))

#
#
## Create grid
#cgrid <- expand.grid(replicate(3, seq(.01, .99, length.out = res), simplify = F))
#cgrid <- cgrid[rowSums(cgrid) < 1,]
#
## keep only cases that could have been observed
#is_obs <- mapply(function(g, rng) g >= rng[1] & g <= rng[2],
#  cgrid, as.data.frame(crngs))
#cgrid <- cgrid[apply(is_obs, 1, all),]

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
  tern_preds <- predict(metamod, newdat)
  
  # Create color
  maxpred <- max(tern_preds, na.rm = T)
  
  colpred <- cols[cut(tern_preds, 
    breaks = seq(-maxpred, maxpred, length.out = length(cols) + 1))
  ]
  attr(colpred, "predictions") <- tern_preds
  colpred
}

tern_values <- TernaryPointValues(tern_pred, direction = 1,
  tern_comp = tern_comp, p = p, ov_mean = ov_mean, metamod = metamod,
  cols = colpal(20), resolution = 100)

# Observed compositions
crngs <- apply(mean_comp[,tern_comp], 2, range)
TPVgrid <- XYToTernary(as.numeric(tern_values[1,]), as.numeric(tern_values[2,]))

is_obs <- mapply(function(g, rng) g >= rng[1] & g <= rng[2],
  as.data.frame(t(TPVgrid)), as.data.frame(crngs))
  
tern_values[3, !apply(is_obs, 1, all)] <- NA

#----- Plot ternary diagram with predictions -----
par(mar = c(5, 4, 4, 7) + .1)
TernaryPlot(alab = spec_labs[tern_comp[1]], blab = spec_labs[tern_comp[2]],
  clab = spec_labs[tern_comp[3]], lab.col = spec_pal[tern_comp],
  grid.col = "white", grid.minor.col = "white")
ColourTernary(tern_values, spectrum = NULL)

# Adding the color scale
coords <- XYToTernary(as.numeric(tern_values[1,]), as.numeric(tern_values[2,]))
real_pred <- attr(tern_pred(coords[1,], coords[2,], coords[3,], 
  tern_comp = tern_comp, p = p, ov_mean = ov_mean, metamod = metamod,
  cols = colpal(20)), "predictions")
maxpred <- max(real_pred)
image.plot(zlim = exp(c(-maxpred, maxpred)), col = colpal(20), 
  legend.only = T)
mtext("RR", side = 3, at = par("usr")[2] + .2, xpd = T, cex = 1.5, line = -1)

dev.print(png, filename = "Results/3_TernaryPredictions.png", 
  units = "in", res = 100)
dev.print(pdf, file = "Results/3_TernaryPredictions.pdf")

#-------------------------------------
#  Diagnostics
#-------------------------------------

# Check interval fo the estimates of variances
# According to the book of Bates, we shouldn't add a random effect for city since we don't have replications.