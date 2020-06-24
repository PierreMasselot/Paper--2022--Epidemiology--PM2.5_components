#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#                using a compositional data approach
#                 with amalgamations
#
#####################################################################

# Amalgamation refers to summing some components. It differs from the balance approach used in papers on regression with compositional data (espcially Hron et al. 2012 and subsequent work). When comparing groups of components, the former compare sums of the components, while the latter compares geometric means. When using amalgamations the subcompositions inside each group don't matter while it does with the balance approach.
# In the first script I wrote (3c), I used balance without really knowing it since papers I read used this approach. However, the issue is that when comparing a component to others, the intra-composition in other matters, while it doesn't with amalgamation. Greenacre (2020) provides a critic of the balance approach.


library(mixmeta)
library(compositions)
library(zCompositions)

load("Data/2_FirstStageResults.RData")

#-------------------------------------
#   Preparing constituent objects
#-------------------------------------

spec_inds <- grep("PM25", colnames(dlist_spec[[1]]))
spec_names <- sapply(strsplit(colnames(dlist_spec[[1]])[spec_inds], "_"),
  "[", 2)
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

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
#   Each components versus all others
#-------------------------------------

# Loop on the constituents to have the isolated contribution of each
res_comp <- matrix(NA, p, 4, 
  dimnames = list(spec_names, c("coef", "lo", "up", "sig")))
for (j in seq_len(p)){
  # Mean summed compositions for each city
  sum_comp <- sapply(dlist_spec, function(x){
    x <- x[,spec_inds]
    # Summing non-focus components
    xs <- cbind(x[,j], rowSums(x[,-j]))
    # Zero value imputation as in Martín-Fernández et al. (2003)
    if (any(xs == 0)) xs <- multRepl(xs, label = 0, dl = rep(1e-5, 2))
    # Compositional mean (pass through the irl transformation)
    mean(acomp(xs))
  })
  sum_comp <- t(sum_comp)
  colnames(sum_comp) <- c(spec_names[j], "Others")
  
  slr <- log(sum_comp[,1] / sum_comp[,2])
  metamod <- mixmeta(coefall ~ slr + indicator, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  
  # We store the SLR coef
  res_comp[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_comp[j, 2] <- res_comp[j, 1] - 1.96 * se_spec
  res_comp[j, 3] <- res_comp[j, 1] + 1.96 * se_spec
  res_comp[j, 4] <- res_comp[j, 2] > 0 | res_comp[j, 3] < 0
}

# Forestplot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_comp[,1], -seq_len(p), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_comp[,2:3]))
segments(res_comp[,2], -seq_len(p), res_comp[, 3], -seq_len(p), 
  col = "darkgrey", lwd = 2)
points(res_comp[,1], -seq_len(p), 
   col = spec_pal, cex = ifelse(res_comp[, 4], 1.5, 1), 
   pch = ifelse(res_comp[, 4], 15, 16))
axis(1)
axis(2, at = -seq_len(p), labels = spec_names, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3e_forestplot_All.png", 
  units = "in", res = 100)
  
#-------------------------------------
#  Aggregated composition
#-------------------------------------

agg_inds <- list(INOR = 1:3, CARBON = 4:5, SS = 6, DUST = 7)
pa <- length(agg_inds)
agg_pal <- c("purple", "darkgreen", 5, 7)

# Loop on the constituents to have the isolated contribution of each
res_agg <- matrix(NA, pa, 4, 
  dimnames = list(names(agg_inds), c("coef", "lo", "up", "sig")))
for (j in seq_len(pa)){
  # Mean summed compositions for each city
  sum_comp <- sapply(dlist_spec, function(x){
    x <- x[,spec_inds]
    # Summing non-focus components
    xs <- cbind(rowSums(x[,agg_inds[[j]], drop = F]), 
      rowSums(x[,-agg_inds[[j]], drop = F]))
    # Zero value imputation as in Martín-Fernández et al. (2003)
    if (any(xs == 0)) xs <- multRepl(xs, label = 0, dl = rep(1e-5, 2))
    # Compositional mean (pass through the irl transformation)
    mean(acomp(xs))
  })
  sum_comp <- t(sum_comp)
  colnames(sum_comp) <- c(names(agg_inds)[j], "Others")
  
  slr <- log(sum_comp[,1] / sum_comp[,2])
  metamod <- mixmeta(coefall ~ slr + indicator, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  
  # We store the SLR coef
  res_agg[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_agg[j, 2] <- res_agg[j, 1] - 1.96 * se_spec
  res_agg[j, 3] <- res_agg[j, 1] + 1.96 * se_spec
  res_agg[j, 4] <- res_agg[j, 2] > 0 | res_agg[j, 3] < 0
}

# Forestplot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_agg[,1], -seq_len(pa), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_agg[,2:3]))
segments(res_agg[,2], -seq_len(pa), res_agg[, 3], -seq_len(pa), 
  col = "darkgrey", lwd = 2)
points(res_agg[,1], -seq_len(pa), 
   col = agg_pal, cex = ifelse(res_agg[, 4], 1.5, 1), 
   pch = ifelse(res_agg[, 4], 15, 16))
axis(1)
axis(2, at = -seq_len(pa), labels = names(agg_inds), 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3e_forestplot_aggregated.png", 
  units = "in", res = 100)
  
#-------------------------------------
#  Fossil-fuel/traffic related subcomposition
#------------------------------------- 
# Sulfate originate from sulfur that mainly comes from fossil fuel combustion (power plants)
# Nitrate is secondary from NOx emission mainly traffic related and gasoline power
# BC is originate from diesel vehicles and traffic more generally
traff_inds <- c(1, 3, 4)

sum_comp <- sapply(dlist_spec, function(x){
  x <- x[,spec_inds]
  # Summing non-focus components
  xs <- cbind(rowSums(x[,traff_inds, drop = F]), 
    rowSums(x[,-traff_inds, drop = F]))
  # Zero value imputation as in Martín-Fernández et al. (2003)
  if (any(xs == 0)) xs <- multRepl(xs, label = 0, dl = rep(1e-5, 2))
  # Compositional mean (pass through the irl transformation)
  mean(acomp(xs))
})
sum_comp <- t(sum_comp)
colnames(sum_comp) <- c("Traffic", "Others")

slr <- log(sum_comp[,1] / sum_comp[,2])
metamod <- mixmeta(coefall ~ slr + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)