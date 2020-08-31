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

# baseline component (actually it doesn't change anything)
ibase <- 1

# Model with ALR as meta-predictors
metamod <- mixmeta(coefall ~ alr(mean_comp, ivar = ibase) + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
  
summary(metamod) # To obtain Q and I2 statistics

# Coefficient retrieving
coef_all <- rep(NA, p)
coef_all[-ibase] <- coef(metamod)[2:7]
coef_all[ibase] <- -sum(coef_all, na.rm = T) # beta7 = - sum(beta1, ..., beta6)

# Standard errors
se_all <- rep(NA, p)
se_all[-ibase] <- sqrt(diag(vcov(metamod))[2:7])
se_all[ibase] <- sqrt(sum(vcov(metamod)[2:7,2:7])) # var(beta7) = sum_j sum_k cov(betaj, betak) (variance of sum of random variables)

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

# Compute predictions 
#preds <- sapply(coef_all, function(b){
#  exp(log(cseq) * b / log(2))
#})

# Compute confience limits
#plo <- sapply(lo_all, function(b){
#  exp(log(cseq) * b / log(2))
#})
#pup <- sapply(up_all, function(b){
#  exp(log(cseq) * b / log(2))
#})

# Prepare objects to store predictions and confidence limits
preds <- plo <- pup <- matrix(NA, length(cseq), p)
for (j in seq_len(p)){
  # Create a compositional grid. The component of interests varies from 0 to 1 (excluded) and the other are taken as the overall mean adjusted for closure while keeping the subcomposition constant
  adj_mean <- sapply(ov_mean[-j] / sum(ov_mean[-j]), "*", 1 - cseq)
  xj <- cbind(cseq, adj_mean)
  
  # Predictions with these data
  bj <- c(coef_all[j], coef_all[-j]) / log(2)
  preds[,j] <- log(xj) %*% bj + coef(metamod)[1]
  
  # Confidence limits
  blj <- c(lo_all[j], coef_all[-j]) / log(2)
  plo[,j] <- log(xj) %*% blj + coef(metamod)[1]
  buj <- c(up_all[j], coef_all[-j]) / log(2)
  pup[,j] <- log(xj) %*% buj + coef(metamod)[1]
}

#--- Plot all prediction together ----
x11()
# Create empty plot
plot(0, 0, col = "white", xlim = c(0,1), ylim = exp(c(min(plo), max(pup))), 
  xlab = "Component proportion", ylab = "RR")
# add polygons for confidence intervals
for (j in seq_len(p)){
  polygon(c(cseq, rev(cseq)), exp(c(plo[,j], rev(pup[,j]))), 
    col = adjustcolor(spec_pal[j], .2), border = NA)
}
# Add predictions
matlines(cseq, exp(preds), lty = 1, lwd = 2, col = spec_pal)
abline(h = 1)
legend("topleft", spec_labs, fill = spec_pal, border = NA, bty = "n", ncol = 2)

dev.print(png, filename = "Results/3_RRpredictions.png", 
  units = "in", res = 200)

#-------------------------------------
#  Ternary plot
#-------------------------------------