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

# baseline component
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

# Forestplot

x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(coef_all, -seq_len(7), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(c(lo_all, up_all)))
abline(v = 0)
segments(lo_all, -seq_len(7), up_all, -seq_len(7), 
  col = "darkgrey", lwd = 2)
points(coef_all, -seq_len(7), cex = ifelse(sig_all, 1.5, 1), 
   pch = ifelse(sig_all, 15, 16), col = spec_pal)
axis(1)
axis(2, at = -seq_len(7), labels = spec_names, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)

dev.print(png, filename = "Paper_Figures/Figure3.png", 
  units = "in", res = 100)

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
  